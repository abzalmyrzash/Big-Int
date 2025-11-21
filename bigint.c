#include "bigint.h"
#include "elog.h"
#include <immintrin.h>
#include <intrin.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <errno.h>
#include <inttypes.h>
#include <limits.h>

int bigint_errno = 0;

typedef BigInt_DataBlock  DataBlock;
typedef BigInt_CapField   CapField;
typedef BigInt_PointField PointField;
typedef BigInt_FormatSpec FormatSpec;

#define CAP_WIDTH   BIGINT_CAP_WIDTH
#define POINT_WIDTH BIGINT_POINT_WIDTH

#define MAX_CAP     BIGINT_MAX_CAP

constexpr DataBlock BLOCK_MAX_VALUE = (DataBlock)~0;

#define HEX_PREFIX "0x"

#define IS_SIGNED_TYPE(TYPE)    !IS_UNSIGNED_TYPE(TYPE)
#define IS_UNSIGNED_TYPE(TYPE)  ((TYPE)-1 > 0)
#define MAX_DEC_INT_DIGITS(TYPE)                \
  ((((sizeof(TYPE) * CHAR_BIT * 1233) >> 12) + 1) \
    + IS_SIGNED_TYPE(TYPE))

#define HEX_DIGITS_PER_BLOCK (sizeof(DataBlock) * 2)
#define MAX_DEC_DIGITS_PER_BLOCK MAX_DEC_INT_DIGITS(DataBlock)
#define DEC_DIGITS_PER_SPACE 3

constexpr int HEX_PREFIX_LEN = sizeof(HEX_PREFIX) - 1;

constexpr size_t BLOCK_WIDTH = sizeof(DataBlock) * CHAR_BIT; // size of DataBlock in bits

static inline uint8_t _addcarry_u8(uint8_t cf, uint8_t a, uint8_t b, uint8_t* out) {
	uint16_t res = (uint16_t)a + b + cf;
	*out = res;
	return res >> 8;
}

static inline uint8_t _subborrow_u8(uint8_t cf, uint8_t a, uint8_t b, uint8_t* out) {
	uint16_t res = (uint16_t)a - b - cf;
	*out = res;
	return res > UINT8_MAX;
}

static inline uint16_t _addcarry_u16(uint16_t cf, uint16_t a, uint16_t b, uint16_t* out) {
	uint32_t res = (uint32_t)a + b + cf;
	*out = res;
	return res >> 16;
}

static inline uint16_t _subborrow_u16(uint16_t cf, uint16_t a, uint16_t b, uint16_t* out) {
	uint32_t res = (uint32_t)a - b - cf;
	*out = res;
	return res > UINT16_MAX;
}

#define CLZ(x) __builtin_clzg(x)
#if BIGINT_BLOCK_WIDTH == 64
	#define ADDCARRY _addcarry_u64
	#define SUBBORROW _subborrow_u64
	typedef __uint128_t BigMult;
	#define FORMAT_HEX     "%" PRIX64
	#define FORMAT_0HEX "%016" PRIX64
	#define FORMAT_hex     "%" PRIx64
	#define FORMAT_0hex "%016" PRIx64

#elif BIGINT_BLOCK_WIDTH == 32
	#define ADDCARRY _addcarry_u32
	#define SUBBORROW _subborrow_u32
	typedef uint64_t BigMult;
	#define FORMAT_HEX     "%" PRIX32
	#define FORMAT_0HEX  "%08" PRIX32
	#define FORMAT_hex     "%" PRIx32
	#define FORMAT_0hex  "%08" PRIx32

#elif BIGINT_BLOCK_WIDTH == 16
	#define ADDCARRY _addcarry_u16
	#define SUBBORROW _subborrow_u16
	typedef uint32_t BigMult;
	#define FORMAT_HEX     "%" PRIX16
	#define FORMAT_0HEX  "%04" PRIX16
	#define FORMAT_hex     "%" PRIx16
	#define FORMAT_0hex  "%04" PRIx16

#elif BIGINT_BLOCK_WIDTH == 8
	#define ADDCARRY _addcarry_u8
	#define SUBBORROW _subborrow_u8
	typedef uint16_t BigMult;
	#define FORMAT_HEX     "%" PRIX8
	#define FORMAT_0HEX  "%02" PRIX8
	#define FORMAT_hex     "%" PRIx8
	#define FORMAT_0hex  "%02" PRIx8
#endif

// returns the low part of a * b, high part in *high
static inline DataBlock block_mul(DataBlock a, DataBlock b, DataBlock* high) {
	BigMult mul = (BigMult)a * b;
	*high = mul >> BLOCK_WIDTH;
	return (DataBlock)mul;
}

struct bigint {
	CapField size;
	CapField cap;
#ifdef BIGINT_SPLIT_POINT_AND_SIGN
	PointField point : POINT_WIDTH - 1;
	PointField sign  : 1;
#else
	PointField point;
	bool sign;
#endif
	DataBlock  data[];
}
#ifdef BIGINT_PACKED
	__attribute__((packed))
#endif
;

constexpr uint8_t DATA_OFFSET = offsetof(struct bigint, data); // offset of data in bytes

static unsigned long long cnt_alloc = 0;
static unsigned long long cnt_realloc = 0;
static unsigned long long cnt_free = 0;

BigInt bigint_bit_set(BigInt z, size_t bitpos);
BigInt bigint_bit_unset(BigInt z, size_t bitpos);
BigInt bigint_bit_toggle(BigInt z, size_t bitpos);

BigInt bigint_alloc(size_t cap) {
	if (cap > MAX_CAP) {
		ELOG_STR("cap exceeds MAX_CAP");
		ELOG("Parameters: new_size = %zu\n", cap);
		return NULL;
	}
	cnt_alloc++;
	BigInt z = malloc(DATA_OFFSET + cap * sizeof(DataBlock));
	if (z == NULL) {
		bigint_errno = BIGINT_ERR_ALLOC_FAIL;
		ELOG_CODE(errno);
		ELOG("Parameters: cap = %zu\n", cap);
		return NULL;
	}
	z->cap = cap;
	z->size = 0;
	z->point = 0;
	z->sign = 0;
	return z;
}

BigInt bigint_zalloc(size_t cap) {
	if (cap > MAX_CAP) {
		ELOG_STR("cap exceeds MAX_CAP");
		ELOG("Parameters: new_size = %zu\n", cap);
		return NULL;
	}
	cnt_alloc++;
	BigInt z = calloc(1, DATA_OFFSET + cap * sizeof(DataBlock));
	if (z == NULL) {
		bigint_errno = BIGINT_ERR_ALLOC_FAIL;
		ELOG_CODE(errno);
		ELOG("Parameters: cap = %zu\n", cap);
		return NULL;
	}
	z->cap = cap;
	return z;
}

BigInt bigint_realloc(BigInt* z_ptr, size_t cap) {
	if (cap > MAX_CAP) {
		ELOG_STR("cap exceeds MAX_CAP");
		ELOG("Parameters: new_size = %zu\n", cap);
		return NULL;
	}
	if (*z_ptr == NULL) return *z_ptr = bigint_alloc(cap);
	cnt_realloc++;
	BigInt z = realloc(*z_ptr, DATA_OFFSET + cap * sizeof(DataBlock));
	if (z == NULL) {
		bigint_errno = BIGINT_ERR_ALLOC_FAIL;
		ELOG_CODE(errno);
		ELOG("Parameters: cap = %zu\n", cap);
		return NULL;
	}
	z->cap = cap;
	return *z_ptr = z;
}

BigInt bigint_rezalloc(BigInt* z_ptr, size_t cap) {
	if (cap > MAX_CAP) {
		ELOG_STR("cap exceeds MAX_CAP");
		ELOG("Parameters: new_size = %zu\n", cap);
		return NULL;
	}
	BigInt z = *z_ptr;
	if (z && z->cap >= cap) {
		memset(z->data, 0, cap * sizeof(DataBlock));
		z->size = 0;
		z->point = 0;
		z->sign = 0;
	}
	else {
		bigint_free(z);
		z = bigint_zalloc(cap);
		if (!z) return NULL;
		*z_ptr = z;
	}
	return z;
}

void bigint_free(BigInt z) {
	if (z) {
		free(z);
		cnt_free++;
	}
}

void bigint_structinfo() {
	printf("sizeof(bigint): %zu\n", sizeof(struct bigint));
	printf("offsetof(bigint, data): %zu\n", offsetof(struct bigint, data));
	printf("sizeof(DataBlock): %zu\n", sizeof(DataBlock));
}

void bigint_memstat() {
	printf("Alloc: %llu\n", cnt_alloc);
	printf("Realloc: %llu\n", cnt_realloc);
	printf("Free: %llu\n", cnt_free);
}

size_t bigint_size(ConstBigInt z) {
	return z->size;
}

size_t bigint_cap(ConstBigInt z) {
	return z->cap;
}

const DataBlock* bigint_data(ConstBigInt z) {
	return z->data;
}

BigInt bigint_resize(BigInt* z_ptr, size_t new_size) {
	if (new_size > MAX_CAP) {
		ELOG_STR("new_size exceeds MAX_CAP");
		ELOG("Parameters: new_size = %zu\n", new_size);
		return NULL;
	}
	BigInt z = *z_ptr;
	if (!z || z->cap < new_size) {
		z = bigint_realloc(z_ptr, new_size);
		if (!z) return NULL;
	}
	/*
	if (new_size == 0) {
		z->point = 0;
		z->sign  = 0;
	}
	*/
	z->size = new_size;
	return z;
}

bool bigint_sign(ConstBigInt z) {
	return z->sign;
}

size_t bigint_point(ConstBigInt z) {
	return z->point;
}

BigInt bigint_abs(ConstBigInt z, BigInt* out_ptr) {
	BigInt out = *out_ptr;
	if (out != z) out = bigint_copy(out_ptr, z);
	out->sign = 0;
	return out;
}

BigInt bigint_neg(ConstBigInt z, BigInt* out_ptr) {
	BigInt out = *out_ptr;
	if (out != z) out = bigint_copy(out_ptr, z);
	out->sign = !z->sign;
	return out;
}

static inline USmallInt ABS(SmallInt v) {
	if (v < 0) return ~(USmallInt)v + 1;
	return (USmallInt)v;
}

BigInt bigint_create_zero() {
	return bigint_zalloc(0);
}

BigInt bigint_set_zero(BigInt* z_ptr) {
	return bigint_resize(z_ptr, 0);
}

BigInt bigint_set_small(BigInt* z_ptr, SmallInt v) {
	if (v == 0) return bigint_set_zero(z_ptr);
	BigInt z = bigint_resize(z_ptr, 1);
	if (!z) return NULL;
	z->data[0] = ABS(v);
	z->point = 0;
	z->sign = (v < 0);
	return z;
}

BigInt bigint_set_usmall(BigInt* z_ptr, USmallInt v) {
	if (v == 0) return bigint_set_zero(z_ptr);
	BigInt z = bigint_resize(z_ptr, 1);
	if (!z) return NULL;
	z->data[0] = v;
	z->point = 0;
	z->sign = 0;
	return z;
}

BigInt bigint_create_small(SmallInt v) {
	if (v == 0) return bigint_create_zero();
	BigInt z = bigint_alloc(1);
	z->data[0] = ABS(v);
	z->size = 1;
	z->point = 0;
	z->sign = (v < 0);
    return z;
}

BigInt bigint_create_usmall(USmallInt v) {
	if (v == 0) return bigint_create_zero();
	BigInt z = bigint_alloc(1);
	z->data[0] = v;
	z->size = 1;
	z->point = 0;
	z->sign = 0;
    return z;
}

BigInt bigint_copy(BigInt* dst_ptr, ConstBigInt src) {
	assert(src);
	assert(dst_ptr);
	BigInt dst = *dst_ptr;
	if (dst == src) return dst;

	if (dst == NULL) dst = bigint_alloc(src->size);
	else if (dst->cap < src->size) {
		bigint_free(dst);
		dst = bigint_alloc(src->size);
	}
	if (dst == NULL) return NULL;
	memcpy(dst->data, src->data, src->size * sizeof(DataBlock));
	dst->point = src->point;
	dst->sign = src->sign;
	dst->size = src->size;
	return *dst_ptr = dst;
}

/*
 * compares the absolute values of a and b
 * 
 * @return value > 0 if |a| > |b|; = 0 if |a| = |b|; < 0 if |a| < |b|
*/
int bigint_ucmp(ConstBigInt a, ConstBigInt b) {
	assert(a);
	assert(b);

	if (a->size > b->size) return 1;
	if (a->size < b->size) return -1;
	if (a->size == 0) return 0;
	size_t i = a->size - 1;
	while (1) {
		if (a->data[i] > b->data[i]) return 1;
		if (a->data[i] < b->data[i]) return -1;
		if (i == 0) break;
		i--;
	}
	return 0;
}

/*
 * compares a and b
 * @return value > 0 if a > b; = 0 if a = b; < 0 if a < b
*/
int bigint_cmp(ConstBigInt a, ConstBigInt b) {
	assert(a);
	assert(b);

	if (a->size == 0 && b->size == 0) return 0;

	if (a->sign && !b->sign) return -1;
	if (!a->sign && b->sign) return 1;
	if (a->sign) return bigint_ucmp(b, a);
	return bigint_ucmp(a, b);
}

USmallInt bigint_usmall(ConstBigInt a) {
	return a->data[0];
}

SmallInt bigint_small(ConstBigInt a) {
	if (a->sign) return -(SmallInt)a->data[0];
	return (SmallInt)a->data[0];
}

int bigint_ucmp_small(ConstBigInt a, USmallInt b) {
	if (a->size == 0) return b;
	if (a->size > 1) return 1;
	return a->data[0] - b;
}

int bigint_cmp_small(ConstBigInt a, SmallInt b) {
	if (a->size == 0) return b;
	if (a->size > 1) return (a->sign) ? -1 : 1;
	return bigint_small(a) - b;
}

// |out| = |a| + |b|
BigInt bigint_uadd(ConstBigInt a, ConstBigInt b, BigInt* out_ptr) {
	assert(a);
	assert(b);
	assert(out_ptr);

	size_t max_size, min_size;
	const DataBlock* max_num;
 
	if (a->size > b->size) {
		max_size = a->size;
		min_size = b->size;
		max_num = a->data;
	}
	else {
		max_size = b->size;
		min_size = a->size;
		max_num = b->data;
	}

	BigInt out;
	if (!(out = bigint_resize(out_ptr, max_size))) return NULL;

	unsigned char carry = 0;
	size_t i = 0;

	for (; i < min_size; i++) {
		carry = ADDCARRY(carry, a->data[i], b->data[i], &out->data[i]);
	}

	for (; i < max_size; i++) {
		carry = ADDCARRY(carry, max_num[i], 0ULL, &out->data[i]);
	}

	if (carry) {
		if (!(out = bigint_resize(out_ptr, max_size + 1))) return NULL;
		out->data[i] = carry;
	}

	return out;
}

// |out| = |a| - |b|,
// assumes |a| >= |b|
BigInt bigint_usub(ConstBigInt a, ConstBigInt b, BigInt* out_ptr) {
	assert(a);
	assert(b);
	assert(out_ptr);
	assert(bigint_ucmp(a, b) >= 0);

	BigInt out;
	if (!(out = bigint_resize(out_ptr, a->size))) return NULL;

	unsigned char borrow = 0;
	size_t i = 0;

	for (; i < b->size; i++) {
		borrow = SUBBORROW(borrow, a->data[i], b->data[i], &out->data[i]);
	}

	for (; i < a->size; i++) {
		borrow = SUBBORROW(borrow, a->data[i], 0ULL, &out->data[i]);
	}

	while (i > 0) {
		if (out->data[i - 1] != 0) break;
		i--;
	}
	out->size = i;

	return out;
}

// out = a + b
BigInt bigint_add(ConstBigInt a, ConstBigInt b, BigInt* out_ptr) {
	assert(a);
	assert(b);
	assert(out_ptr);

	if (a->sign == b->sign) {
		bigint_uadd(a, b, out_ptr);
		(*out_ptr)->sign = a->sign;
	}

	else {
		int cmp = bigint_ucmp(a, b);
		if (cmp > 0) { // |a| > |b|
			bigint_usub(a, b, out_ptr);
			(*out_ptr)->sign = a->sign;
		}
		else if (cmp < 0) { // |a| < |b|
			bigint_usub(b, a, out_ptr);
			(*out_ptr)->sign = b->sign;
		}
		else { // |a| = |b|
			if(*out_ptr == NULL) *out_ptr = bigint_create_zero();
			else bigint_set_zero(out_ptr);
		}
	}

	return *out_ptr;
}

// out = a - b,
BigInt bigint_sub(ConstBigInt a, ConstBigInt b, BigInt* out_ptr) {
	assert(a);
	assert(b);
	assert(out_ptr);

	if (a->sign != b->sign) {
		bigint_uadd(a, b, out_ptr);
		(*out_ptr)->sign = a->sign;
	}

	else {
		int cmp = bigint_ucmp(a, b);
		if (cmp > 0) { // |a| > |b|
			bigint_usub(a, b, out_ptr);
			(*out_ptr)->sign = a->sign;
		}
		else if (cmp < 0) { // |a| < |b|
			bigint_usub(b, a, out_ptr);
			(*out_ptr)->sign = !a->sign;
		}
		else { // |a| = |b|
			if(*out_ptr == NULL) *out_ptr = bigint_create_zero();
			else bigint_set_zero(out_ptr);
		}
	}

	return *out_ptr;
}

// |out| = |a| * |b|
BigInt bigint_umul(ConstBigInt a, ConstBigInt b, BigInt* out_ptr) {
	assert(a);
	assert(b);
	assert(out_ptr);

	size_t size = a->size + b->size;

	BigInt out = *out_ptr;

	BigInt to_be_freed = NULL;
	if (out == a) {
		BigInt tmp = NULL;
		bigint_copy(&tmp, a);
		if (!tmp) return NULL;
		a = tmp;
		to_be_freed = (BigInt)a;
	}
	else if (out == b) {
		BigInt tmp = NULL;
		bigint_copy(&tmp, b);
		if (!tmp) return NULL;
		b = tmp;
		to_be_freed = (BigInt)b;
	}
	if (!bigint_rezalloc(&out, size)) goto error;

	DataBlock* arr = out->data;

	for (size_t i = 0; i < a->size; i++) {
		DataBlock d1 = a->data[i];
		for (size_t j = 0; j < b->size; j++) {
			DataBlock d2 = b->data[j];
			DataBlock high;
			DataBlock low = block_mul(d1, d2, &high);
			unsigned char carry;
			carry = ADDCARRY(0, arr[i + j], low, &arr[i + j]);
			carry = ADDCARRY(carry, arr[i + j + 1], high, &arr[i + j + 1]);
			size_t idx = i + j + 2;
			while (carry) {
				carry = ADDCARRY(carry, arr[idx], 0, &arr[idx]);
				idx++;
			}
		}
	}

	while (size > 0) {
		if (arr[size - 1] != 0) break;
		size--;
	}
	out->size = size;

	bigint_free(to_be_freed);
	return *out_ptr = out;
error:
	bigint_free(to_be_freed);
	return NULL;
}

// out = a * b
BigInt bigint_mul(ConstBigInt a, ConstBigInt b, BigInt* out_ptr) {
	assert(a);
	assert(b);
	assert(out_ptr);

	bool out_sign = a->sign ^ b->sign;
	bigint_umul(a, b, out_ptr);
	(*out_ptr)->sign = out_sign;
	return *out_ptr;
}

// |out| = |a| / |b|,
// |rem| = |a| % |b|
BigInt bigint_udiv(ConstBigInt a, ConstBigInt b, BigInt* out_ptr, BigInt* rem_ptr) {
	assert(a);
	assert(b);
	assert(out_ptr);
	if (b->size == 0) {
		bigint_errno = BIGINT_ERR_DIV_BY_ZERO;
		ELOG_STR("DIVISION BY ZERO");
		return NULL;
	}

	BigInt out = *out_ptr;
	BigInt rem = NULL;
	BigInt c = NULL;

	bool a_copy = false;
	bool b_copy = false;

	if (out == a) {
		BigInt tmp = NULL;
		bigint_copy(&tmp, a);
		if (!tmp) goto error;
		a = tmp;
		a_copy = true;
	}

	else if (out == b) {
		BigInt tmp = NULL;
		bigint_copy(&tmp, b);
		if (!tmp) goto error;
		b = tmp;
		b_copy = true;
	}

	if (rem_ptr) {
		rem = *rem_ptr;
		if (rem == a) {
			assert(!a_copy);
			BigInt tmp = NULL;
			bigint_copy(&tmp, a);
			if (!tmp) goto error;
			a = tmp;
			a_copy = true;
		}

		else if (rem == b) {
			assert(!b_copy);
			BigInt tmp = NULL;
			bigint_copy(&tmp, b);
			if (!tmp) goto error;
			b = tmp;
			b_copy = true;
		}
		assert(rem == NULL || rem != out);
	}

	if (bigint_ucmp_small(b, 1) == 0) {
		bigint_copy(&out, a);
		if (!out) goto error;
		if (rem_ptr) bigint_set_zero(&rem);
		goto ret;
	}

	int cmp = bigint_ucmp(a, b);

	if (cmp < 0) {
		if (rem_ptr) {
			bigint_copy(&rem, a);
			if (!rem) goto error;
		}
		bigint_set_zero(&out);
		if (!out) goto error;
	}

	else if (cmp == 0) {
		if (rem_ptr) {
			bigint_set_zero(&rem);
			if (!rem) goto error;
		}
		bigint_set_usmall(&out, 1);
		if (!out) goto error;
	}
	
	else {
		size_t size_diff = a->size - b->size;
		size_t bitpos;
		bigint_copy(&rem, a);
		if (!rem) goto error;
		if (!bigint_lshift_digits(b, size_diff, &c)) goto error;

		if (bigint_ucmp(rem, c) < 0) {
			if (!bigint_rezalloc(&out, size_diff)) goto error;
			out->size = size_diff;
			bigint_rshift_digits(c, 1, &c);
			bitpos = (size_diff - 1) * BLOCK_WIDTH;
		} else {
			if (!bigint_rezalloc(&out, size_diff + 1)) goto error;
			out->size = size_diff + 1;
			bitpos = size_diff * BLOCK_WIDTH;
		}

		do {
			if (!bigint_lshift(c, 1, &c)) goto error;
			bitpos++;
		} while (bigint_ucmp(rem, c) >= 0);

		bigint_rshift(c, 1, &c);
		bigint_usub(rem, c, &rem);
		bigint_bit_set(out, --bitpos);

		while(bigint_ucmp(rem, b) >= 0) {
			while (bigint_ucmp(rem, c) < 0) {
				bigint_rshift(c, 1, &c);
				assert(bitpos > 0);
				bitpos--;
			}
			bigint_usub(rem, c, &rem);
			bigint_bit_set(out, bitpos);
		}
	}

ret:
	/*
	BigInt check = NULL;
	bigint_umul(out, b, &check);
	assert(check);
	bigint_uadd(check, rem, &check);
	assert(check);
	assert(bigint_ucmp(check, a) == 0 || bigint_printf("\n%px\n\n", check));
	free(check);
	*/
	bigint_free(c);
	if (a_copy) bigint_free((BigInt)a);
	if (b_copy) bigint_free((BigInt)b);
	if (rem_ptr) {
		*rem_ptr = rem;
	} else {
		bigint_free(rem);
	}
	return *out_ptr = out;

error:
	bigint_free(c);
	if (a_copy) bigint_free((BigInt)a);
	if (b_copy) bigint_free((BigInt)b);
	if (rem_ptr) {
		*rem_ptr = rem;
	} else {
		bigint_free(rem);
	}
	return NULL;
}

// out = a / b,
// rem = a % b
BigInt bigint_div(ConstBigInt a, ConstBigInt b, BigInt* out_ptr, BigInt* rem_ptr) {
	assert(a);
	assert(b);
	assert(out_ptr);
	bool out_sign = a->sign ^ b->sign;
	bool rem_sign = a->sign;
	if (!bigint_udiv(a, b, out_ptr, rem_ptr)) return NULL;
	(*out_ptr)->sign = out_sign;
	if (rem_ptr) (*rem_ptr)->sign = rem_sign;
	return *out_ptr;
}

BigInt bigint_lshift_digits(ConstBigInt z, size_t shift, BigInt* out_ptr) {
	assert(z);
	assert(out_ptr);
	if (z->size == 0) return bigint_set_zero(out_ptr);
	if (shift == 0) return bigint_copy(out_ptr, z);

	size_t out_size = z->size + shift;
	if (!bigint_resize(out_ptr, out_size)) return NULL;
	memset((*out_ptr)->data, 0, shift * sizeof(DataBlock));
	memmove((*out_ptr)->data + shift, z->data, z->size * sizeof(DataBlock));
	return *out_ptr;
}

BigInt bigint_rshift_digits(ConstBigInt z, size_t shift, BigInt* out_ptr) {
	assert(z);
	assert(out_ptr);
	if (z->size <= shift) return bigint_set_zero(out_ptr);
	if (shift == 0) return bigint_copy(out_ptr, z);

	size_t out_size = z->size - shift;
	if (!bigint_resize(out_ptr, out_size)) return NULL;
	memmove((*out_ptr)->data, z->data + shift, out_size * sizeof(DataBlock));
	return *out_ptr;
}

BigInt bigint_lshift(ConstBigInt z, size_t shift, BigInt* out_ptr) {
	assert(z);
	assert(out_ptr);
	if (z->size == 0) return bigint_set_zero(out_ptr);

	const size_t shift_digits = shift / BLOCK_WIDTH;
	const size_t shift_bits = shift % BLOCK_WIDTH;
	if (shift_bits == 0) return bigint_lshift_digits(z, shift_digits, out_ptr);

	const DataBlock LOW_BITS = (1ULL << (BLOCK_WIDTH - shift_bits)) - 1;
	const DataBlock HIGH_BITS = ~0ULL - LOW_BITS;

	size_t out_size = z->size + shift_digits;

	size_t i = z->size - 1;
	DataBlock hi, lo;

	hi = z->data[i] & HIGH_BITS;
	lo = z->data[i] & LOW_BITS;
	BigInt out;
	if (hi) {
		out = bigint_resize(out_ptr, ++out_size);
		if (!out) return NULL;
		out->data[i + shift_digits + 1] = hi >> (BLOCK_WIDTH - shift_bits);
	} else {
		out = bigint_resize(out_ptr, out_size);
		if (!out) return NULL;
	}
	out->data[i + shift_digits] = lo << shift_bits;

	if (i--) while (1) {
		hi = z->data[i] & HIGH_BITS;
		lo = z->data[i] & LOW_BITS;
		out->data[i + shift_digits] = lo << shift_bits;
		out->data[i + shift_digits + 1] |= hi >> (BLOCK_WIDTH - shift_bits);
		if (i == 0) break;
		i--;
	}

	memset(out->data, 0, shift_digits * sizeof(DataBlock));
	return out;
}

BigInt bigint_rshift(ConstBigInt z, size_t shift, BigInt* out_ptr) {
	assert(z);
	assert(out_ptr);

	size_t shift_digits = shift / BLOCK_WIDTH;
	size_t shift_bits = shift % BLOCK_WIDTH;
	if (shift_bits == 0) return bigint_rshift_digits(z, shift_digits, out_ptr);

	size_t z_size = z->size;
	if (shift_digits >= z_size) {
		return bigint_set_zero(out_ptr);
	}

	size_t out_size = z_size - shift_digits;
	if (shift_bits >= BLOCK_WIDTH - CLZ((DataBlock)z->data[z_size - 1])) {
		if (out_size <= 1) return bigint_set_zero(out_ptr);
		out_size--;
	}
	BigInt out = bigint_resize(out_ptr, out_size);
	if (!out) return NULL;

	const DataBlock LOW_BITS = (1ULL << shift_bits) - 1;
	const DataBlock HIGH_BITS = ~0ULL - LOW_BITS;

	size_t i;
	DataBlock hi, lo;
	for (i = 0; i < out_size - 1; i++) {
		hi = z->data[i + shift_digits] & HIGH_BITS;
		lo = z->data[i + shift_digits + 1] & LOW_BITS;
		out->data[i] = (hi >> shift_bits) | lo << (BLOCK_WIDTH - shift_bits);
	}

	hi = z->data[i + shift_digits] & HIGH_BITS;
	if (i + shift_digits + 1 < z_size)
		lo = z->data[i + shift_digits + 1] & LOW_BITS;
	else lo = 0;
	out->data[i] = (hi >> shift_bits) | lo << (BLOCK_WIDTH - shift_bits);

	return out;
}

BigInt bigint_bit_set(BigInt z, size_t bitpos) {
	assert(z);
	z->data[bitpos / BLOCK_WIDTH] |= (1ULL << (bitpos % BLOCK_WIDTH));
	return z;
}

BigInt bigint_bit_unset(BigInt z, size_t bitpos) {
	assert(z);
	z->data[bitpos / BLOCK_WIDTH] &= ~(1ULL << (bitpos % BLOCK_WIDTH));
	return z;
}

BigInt bigint_bit_toggle(BigInt z, size_t bitpos) {
	assert(z);
	z->data[bitpos / BLOCK_WIDTH] ^= (1ULL << (bitpos % BLOCK_WIDTH));
	return z;
}

int bigint_fwrite(FILE* file, ConstBigInt z, bool is_signed) {
	assert(z);
	uint64_t size = z->size;
	if (fwrite(&size, sizeof(size), 1, file) != 1) {
		return errno;
	}
	uint64_t point = z->point;
	if (fwrite(&point, sizeof(point), 1, file) != 1) {
		return errno;
	}
	bool sign = z->sign;
	if (is_signed && fwrite(&sign, sizeof(sign), 1, file) != 1) {
		return errno;
	}
	if (fwrite(&z->data, sizeof(*z->data), z->size, file) != z->size) {
		return errno;
	}
	return 0;
}

int bigint_fread(FILE* file, BigInt* z_ptr, bool is_signed) {
	assert(z_ptr);
	uint64_t size;
	if (fread(&size, sizeof(size), 1, file) != 1) return errno;
	BigInt z = bigint_resize(z_ptr, size);
	if (!z) return errno;

	uint64_t point;
	if (fread(&point, sizeof(point), 1, file) != 1) return errno;
	z->point = point;

	bool sign;
	if (is_signed && fread(&sign, sizeof(sign), 1, file) != 1) return errno;
	z->sign = sign;

	if (fread(&z->data, sizeof(*z->data), size, file) != size) return errno;
	return 0;
}

FormatSpec bigint_initBIFS() {
	return (FormatSpec) { 0 };
}

int _bigint_fprintf(FILE* file, const char* const format, va_list args) {
	char c = *format;
	int i = 0;
	FormatSpec bifs = bigint_initBIFS();
	bool processing_bifs = false;
	int charsWritten = 0;

	while (c != '\0') {
		if (processing_bifs) {
			switch(c) {
				case 'u': case 'U':
					bifs.is_unsigned = true;
					break;

				case '+':
					bifs.add_plus_sign = true;
					break;

				case 'p': case 'P':
					bifs.add_prefix = true;
					break;

				case '0':
					bifs.leading_zeroes = true;
					break;

				case '_':
					bifs.add_spaces = true;
					break;

				case 'x': case 'X':
					bifs.base = 16;
					goto read_big_int;

				case 'd': case 'D':
					bifs.base = 10;

				read_big_int:
					ConstBigInt z = va_arg(args, ConstBigInt);
					if (!z) {
						ELOG_CODE_STR(errno, "NULL argument");
						return -1;
					}

					bifs.uppercase = (c == 'X');
					char* str = bigint_str(z, bifs);
					if (!str) {
						ELOG_CODE_STR(errno, "failed to allocate str");
						return -1;
					}

					int res = fprintf(file, "%s", str);
					if (res < 0) {
						ELOG_CODE_STR(errno, "fprintf failed");
						return -1;
					}

					charsWritten += res;
					free(str);
					processing_bifs = false;
					bifs = bigint_initBIFS();
					break;

				case '%':
					int err = putc('%', file);
					if (res < 0) {
						ELOG_CODE_STR(errno, "putc failed");
						return res;
					}
					processing_bifs = false;
					bifs = bigint_initBIFS();
					break;

				default:
					fprintf(stderr, "ERROR in _bigint_fprint(%s): Invalid format specifier.\n",
							format);
					return -1;
					break;
			}
		}

		else if (c == '%') {
			processing_bifs = true;
		}
		else {
			int res = putc(c, file);
			if (res < 0) {
				fprintf(stderr, "ERROR in %s:%d:%s(): putc failed with code %d.\n",
						__FILE__, __LINE__, __FUNCTION__, res);
				return res;
			}
			charsWritten++;
		}

		c = format[++i];
	}

	return charsWritten;
}

int bigint_fprintf(FILE* file, const char* const format, ...) {
	va_list args;
	va_start(args, format);
	int res = _bigint_fprintf(file, format, args);
	va_end(args);
	return res;
}

int bigint_printf(const char* const format, ...) {
	va_list args;
	va_start(args, format);
	int res = _bigint_fprintf(stdout, format, args);
	va_end(args);
	return res;
}

char* bigint_str(ConstBigInt z, FormatSpec bifs) {
	if (bifs.base == 16) return bigint_hex_str(z, bifs);
	if (bifs.base == 10) return bigint_dec_str(z, bifs);
	assert(false && "Invalid base provided!");
	return NULL;
}

char* bigint_hex_str(ConstBigInt z, FormatSpec bifs) {
	assert(z);
	bool add_minus_sign = z->sign && !bifs.is_unsigned; 
	bool add_sign = add_minus_sign || bifs.add_plus_sign;

	const size_t max_digits_and_spaces = 
		(z->size > 0) ? (z->size * (HEX_DIGITS_PER_BLOCK + bifs.add_spaces)) : 1;
	char* str = malloc(max_digits_and_spaces
			+ add_sign + bifs.add_prefix * HEX_PREFIX_LEN + 1);
	if (!str) return NULL;

	const char* const f_0hex = bifs.uppercase ? FORMAT_0HEX : FORMAT_0hex;
	const char* const f_hex = bifs.uppercase ? FORMAT_HEX : FORMAT_hex;

	size_t offset = 0;
	size_t res;

	if (add_minus_sign) {
		str[offset++] = '-';
	} else if (bifs.add_plus_sign) {
		str[offset++] = '+';
	}

	if (bifs.add_prefix) {
		offset += sprintf(str + offset, HEX_PREFIX);
	}

	if (z->size == 0) {
		str[offset++] = '0';
		str[offset] = '\0';
		return str;
	}

	size_t i = z->size - 1;
	if (bifs.leading_zeroes) offset += sprintf(str + offset, f_0hex, z->data[i]);
	else offset += sprintf(str + offset, f_hex, z->data[i]);

	while (i-- > 0) {
		if (bifs.add_spaces) offset += sprintf(str + offset, " ");
		offset += sprintf(str + offset, f_0hex, z->data[i]);
	}

	return str;
}

char* bigint_dec_str(ConstBigInt z, FormatSpec bifs) {
	assert(z);
	const bool add_minus_sign = z->sign && !bifs.is_unsigned; 
	const bool add_sign = add_minus_sign || bifs.add_plus_sign;

	const size_t max_digits =
		(z->size > 0) ? z->size * MAX_DEC_DIGITS_PER_BLOCK : 1;
	const size_t max_spaces = (max_digits / DEC_DIGITS_PER_SPACE) * bifs.add_spaces;
	char* const str = malloc(max_digits + max_spaces + add_sign + 1);
	if (!str) return NULL;

	size_t offset = 0;

	if (add_minus_sign) {
		str[offset++] = '-';
	} else if (bifs.add_plus_sign) {
		str[offset++] = '+';
	}

	if (z->size == 0) {
		str[offset++] = '0';
		str[offset] = '\0';
		return str;
	}

	char* const digits = str + offset;
	size_t cnt = 0;

	BigInt a = NULL;
	if (!bigint_copy(&a, z)) return NULL;

	BigInt ten = bigint_create_usmall(10);
	if (!ten) return NULL;
	BigInt rem = NULL;

	while (a->size > 0) {
		bigint_udiv(a, ten, &a, &rem);
		assert(cnt < max_digits + max_spaces);
		digits[cnt++] = bigint_usmall(rem) + '0';
		if (bifs.add_spaces &&
				(cnt + 1) % (DEC_DIGITS_PER_SPACE + 1) == 0) {
			assert(cnt < max_digits + max_spaces);
			digits[cnt++] = ' ';
		}
	}
	if (digits[cnt - 1] == ' ') cnt--;

	bigint_free(a);
	bigint_free(ten);
	bigint_free(rem);

	for (size_t i = 0; i < cnt / 2; i++) {
		char tmp = digits[i];
		digits[i] = digits[cnt - i - 1];
		digits[cnt - i - 1] = tmp;
	}
	digits[cnt] = '\0';

	return str;
}

BigInt bigint_scan_hex(const char* str, size_t str_len, BigInt* out_ptr);
BigInt bigint_scan_dec(const char* str, size_t str_len, BigInt* out_ptr);

BigInt bigint_scan(const char* str, BigInt* out_ptr) {
	size_t str_len = strlen(str);
	size_t offset = 0;
	bool sign = false;
	BigInt out;

	char signch = str[offset];
	if (signch == '-') {
		sign = true;
		offset++;
	} else if (signch == '+') {
		offset++;
	}

	bool hex = true;
	for (int i = 0; i < HEX_PREFIX_LEN; i++) {
		if (HEX_PREFIX[i] != str[offset + i]) {
			hex = false;
			break;
		}
	}
	if (hex) {
		offset += HEX_PREFIX_LEN;
		out = bigint_scan_hex(str + offset, str_len - offset, out_ptr);
	} else {
		out = bigint_scan_dec(str + offset, str_len - offset, out_ptr);
	}

	if (!out) return NULL;
	out->sign = sign;
	return out;
}

BigInt bigint_scan_hex(const char* str, size_t str_len, BigInt* out_ptr) {
	size_t cap = str_len / HEX_DIGITS_PER_BLOCK + 1;
	BigInt out = bigint_rezalloc(out_ptr, cap);
	if (!out) return NULL;
	size_t shift = 0;

	size_t i = str_len;
	while (i > 0) {
		i--;
		char c = str[i];
		int value;
		if (c >= '0' && c <= '9') {
			value = c - '0';
		}
		else if (c >= 'A' && c <= 'F') {
			value = c - 'A' + 10;
		}
		else if (c >= 'a' && c <= 'f') {
			value = c - 'a' + 10;
		} else {
			continue;
		}

		if (value) {
			out->size = shift / BLOCK_WIDTH + 1;
			if (value & 1)
				bigint_bit_set(out, shift);
			if (value & 2)
				bigint_bit_set(out, shift + 1);
			if (value & 4)
				bigint_bit_set(out, shift + 2);
			if (value & 8)
				bigint_bit_set(out, shift + 3);
		}
		shift += 4;
	}

	return out;
}

BigInt bigint_scan_dec(const char* str, size_t str_len, BigInt* out_ptr) {
	size_t cap = str_len / (MAX_DEC_DIGITS_PER_BLOCK - 1) + 1;
	BigInt out = bigint_rezalloc(out_ptr, cap);
	if (!out) return NULL;

	BigInt ten = bigint_create_usmall(10);
	if (!ten) return NULL;

	BigInt pow10 = bigint_alloc(cap);
	if (!pow10) return NULL;
	bigint_set_usmall(&pow10, 1);

	BigInt digit = bigint_alloc(cap);
	if (!digit) return NULL;

	size_t i = str_len;
	while (i > 0) {
		i--;
		char c = str[i];
		int value;
		if (c >= '0' && c <= '9') {
			value = c - '0';
		}
		else {
			continue;
		}

		bigint_set_usmall(&digit, value);
		bigint_umul(digit, pow10, &digit);
		bigint_uadd(out, digit, &out);

		bigint_umul(pow10, ten, &pow10);
	}
	bigint_free(pow10);
	bigint_free(ten);
	bigint_free(digit);

	return out;
}
