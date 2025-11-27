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

// Meant to be used with an array successively. Assigns no leading zeros.
// Doesn't assign the result to out if it's zero, but increments cnt_zeros.
// If it is non-zero, assigns result to out, sets cnt_zeros blocks
// that came before out, to 0, and sets cnt_zeros back to 0.
static inline uint8_t subborrow_nlz(uint8_t cf, DataBlock a, DataBlock b,
		DataBlock* out, size_t* cnt_zeros) {	
	DataBlock block;
	uint8_t borrow = SUBBORROW(cf, a, b, &block);
	if (block == 0) {
		(*cnt_zeros)++;
	}
	else {
		*out = block;
		if (*cnt_zeros > 0) {
			memset(out - *cnt_zeros, 0, *cnt_zeros * sizeof(DataBlock));
			*cnt_zeros = 0;
		}
	}
	return borrow;
}

// size:  length of number in datablocks
// cap:   number of datablocks allocated
// point: position of the fractional point
// sign:  0 means positive, 1 means negative
struct bigint {
	CapField   size;
	CapField   cap;
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

typedef struct {
	const DataBlock* data;
	CapField size;
} Slice;

typedef struct {
	void* buf;
	size_t bufsize;
} Context;

static thread_local Context ctx;

void bigint_init() {
	ctx.buf = NULL;
	ctx.bufsize = 0;
}

void bigint_finish() {
	free(ctx.buf);
}

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
	BigInt z = malloc(DATA_OFFSET + cap * sizeof(DataBlock));
	if (z == NULL) {
		bigint_errno = BIGINT_ERR_ALLOC_FAIL;
		ELOG_CODE(errno);
		ELOG("Parameters: cap = %zu\n", cap);
		return NULL;
	}
	cnt_alloc++;
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
	BigInt z = calloc(1, DATA_OFFSET + cap * sizeof(DataBlock));
	if (z == NULL) {
		bigint_errno = BIGINT_ERR_ALLOC_FAIL;
		ELOG_CODE(errno);
		ELOG("Parameters: cap = %zu\n", cap);
		return NULL;
	}
	cnt_alloc++;
	z->cap = cap;
	return z;
}

BigInt bigint_realloc(BigInt* z_ptr, size_t cap, bool keep_data) {
	if (cap > MAX_CAP) {
		ELOG_STR("cap exceeds MAX_CAP");
		ELOG("Parameters: new_size = %zu\n", cap);
		return NULL;
	}
	BigInt z = *z_ptr;
	if (!z) return *z_ptr = bigint_alloc(cap);
	if (keep_data || cap < z->cap) {
		z = realloc(z, DATA_OFFSET + cap * sizeof(DataBlock));
		if (z) cnt_realloc++;
	} else {
		bigint_free(z);
		z = bigint_alloc(cap);
	}
	if (z == NULL) {
		bigint_errno = BIGINT_ERR_ALLOC_FAIL;
		ELOG_CODE(errno);
		ELOG("Parameters: cap = %zu\n", cap);
		return NULL;
	}
	z->cap = cap;
	return *z_ptr = z;
}

BigInt bigint_realloc_if_small(BigInt* z_ptr, size_t cap, bool keep_data) {
	if (!(*z_ptr) || (*z_ptr)->cap < cap) {
		return bigint_realloc(z_ptr, cap, keep_data);
	}
	return *z_ptr;
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

static BigInt bigint_resize(BigInt* z_ptr, size_t new_size) {
	BigInt z = *z_ptr;
	z = bigint_realloc_if_small(z_ptr, new_size, true);
	if (!z) return NULL;
	z->size = new_size;
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

static inline size_t bitsize(ConstBigInt z) {
	return z->size * BLOCK_WIDTH - CLZ(z->data[z->size - 1]);
}

size_t bigint_cap(ConstBigInt z) {
	return z->cap;
}

const DataBlock* bigint_data(ConstBigInt z) {
	return z->data;
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
	BigInt z = *z_ptr;
	if (!z) return *z_ptr = bigint_zalloc(0);
	z->size = 0;
	return z;
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

// standard comparison function
static int cmp(Slice a, Slice b) {
	if (a.size > b.size) return 1;
	if (a.size < b.size) return -1;
	if (a.size == 0) return 0;
	size_t i = a.size - 1;
	while (1) {
		if (a.data[i] > b.data[i]) return 1;
		if (a.data[i] < b.data[i]) return -1;
		if (i == 0) return 0;
		i--;
	}
}

// Same as cmp, but you pass pointer to diff_point (can't be NULL).
// For a and b that have the same size,
// diff_point is the point at which a and b differ
// (counting from least significant block).
// For example (let's assume the digits here are blocks):
// a = 12345, b = 12356. Their first three digits are the same,
// leaving the last two digits (45 in a and 56 in b) that are different.
// Therefore diff_point = 2.
// If one number's size is larger than the other, diff_point = the larger size.
// If both numbers are the same, diff_point = 0.
static int cmpd(Slice a, Slice b, size_t* diff_point) {
	assert(diff_point);
	if (a.size > b.size) {
		*diff_point = a.size;
		return 1;
	}
	if (a.size < b.size) {
		*diff_point = b.size;
		return -1;
	}
	if (a.size == 0) {
		*diff_point = 0;
		return 0;
	}
	size_t i = a.size - 1;
	while (1) {
		if (a.data[i] > b.data[i]) {
			*diff_point = i + 1;
			return 1;
		}
		if (a.data[i] < b.data[i]) {
			*diff_point = i + 1;
			return -1;
		}
		if (i == 0) {
			*diff_point = 0;
			return 0;
		}
		i--;
	}
}

// compare (a << lshift) with b
// NOTE: lshift by number of blocks not bits
static int cmp_lshifted(Slice a, size_t lshift, Slice b) {
	if (a.size == 0) {
		if (b.size > 0) return -1;
		return 0;
	}
	if (a.size + lshift > b.size) return 1;
	if (a.size + lshift < b.size) return -1;
	size_t i = b.size - 1;
	while (1) {
		if (a.data[i - lshift] > b.data[i]) return 1;
		if (a.data[i - lshift] < b.data[i]) return -1;
		if (i == lshift) break;
		i--;
	}
	while (1) {
		if (i == 0) return 0;
		i--;
		if (b.data[i] > 0) return -1;
	}
}

// compare (a << lshift) with b
// see cmpd() for info about diff_point
// NOTE: lshift by number of blocks not bits
static int cmpd_lshifted(Slice a, size_t lshift, Slice b, size_t* diff_point) {
	if (a.size == 0) {
		if (b.size > 0) {
			*diff_point = b.size;
			return -1;
		}
		*diff_point = 0;
		return 0;
	}
	if (a.size + lshift > b.size) {
		*diff_point = a.size + lshift;
		return 1;
	}
	if (a.size + lshift < b.size) {
		*diff_point = b.size;
		return -1;
	}
	size_t i = b.size - 1;
	while (1) {
		if (a.data[i - lshift] > b.data[i]) {
			*diff_point = i + 1;
			return 1;
		}
		if (a.data[i - lshift] < b.data[i]) {
			*diff_point = i + 1;
			return -1;
		}
		if (i == lshift) {
			break;
		}
		i--;
	}
	while (1) {
		if (i == 0) {
			*diff_point = 0;
			return 0;
		}
		i--;
		if (b.data[i] > 0) {
			*diff_point = i + 1;
			return -1;
		}
	}
}

/*
 * compares the absolute values of a and b
 * 
 * @return value > 0 if |a| > |b|; = 0 if |a| = |b|; < 0 if |a| < |b|
*/
int bigint_ucmp(ConstBigInt a, ConstBigInt b) {
	assert(a);
	assert(b);
	const Slice A = { .data = a->data, .size = a->size };
	const Slice B = { .data = b->data, .size = b->size };
	return cmp(A, B);
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
	if (a->size == 0) return 0;
	return a->data[0];
}

SmallInt bigint_small(ConstBigInt a) {
	if (a->size == 0) return 0;
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

static void print_slice(Slice a) {
	if (a.size == 0) {
		printf("0");
		return;
	}
	printf(FORMAT_HEX, a.data[a.size - 1]);
	a.size--;
	while(a.size > 0) {
		printf(" " FORMAT_0HEX, a.data[a.size - 1]);
		a.size--;
	}
}

static inline size_t size_without_zeros(const DataBlock* a, size_t size) {
	while (size > 0) {
		if (a[size - 1] != 0) break;
		size--;
	}
	return size;
}

static inline Slice split(Slice* slice, size_t size) {
	Slice second;
	if (slice->size > size) {
		second.data = slice->data + size;
		second.size = size_without_zeros(second.data, slice->size - size);
		slice->size = size_without_zeros(slice->data, size);
	} else {
		second.data = NULL;
		second.size = 0;
	}
	return second;
}

// a.size must be >= b.size
static size_t add_a(Slice a, Slice b, DataBlock* out, bool assign_carry) {
	assert(a.size >= b.size);
	assert(a.size == size_without_zeros(a.data, a.size));
	assert(b.size == size_without_zeros(b.data, b.size));
	int carry = 0;
	size_t i = 0;
	for (; i < b.size; i++) {
		carry = ADDCARRY(carry, a.data[i], b.data[i], &out[i]);
	}
	for (; i < a.size; i++) {
		carry = ADDCARRY(carry, a.data[i], 0, &out[i]);
	}
	if (assign_carry && carry) out[a.size] = carry;
	return a.size + carry;
}

static size_t add(Slice a, Slice b, DataBlock* out, bool assign_carry) {
	if (a.size >= b.size) {
		return add_a(a, b, out, assign_carry);
	}
	else {
		return add_a(b, a, out, assign_carry);
	}
}

// out = (a << lshift) + b;
// NOTE: lshift by number of blocks not bits
static size_t add_lshifted(Slice a, size_t lshift, Slice b, DataBlock* out, bool assign_carry) {
	if (a.size == 0) {
		memmove(out, b.data, b.size * sizeof(DataBlock));
		return b.size;
	}
	if (b.size <= lshift) {
		memmove(out, b.data, b.size * sizeof(DataBlock));
		memset(out + b.size, 0, (lshift - b.size) * sizeof(DataBlock));
		memmove(out + lshift, a.data, a.size * sizeof(DataBlock));
		return a.size + lshift;
	}
	memmove(out, b.data, lshift * sizeof(DataBlock));
	int carry = 0;
	size_t i = lshift;
	if (a.size + lshift > b.size) {
		for (; i < b.size; i++) {
			carry = ADDCARRY(carry, a.data[i - lshift], b.data[i], &out[i]);
		}
		for (; i < a.size + lshift; i++) {
			carry = ADDCARRY(carry, a.data[i - lshift], 0, &out[i]);
		}
		if (assign_carry && carry) out[a.size + lshift] = carry;
		return a.size + lshift + carry;
	} else {
		for (; i < a.size + lshift; i++) {
			carry = ADDCARRY(carry, a.data[i - lshift], b.data[i], &out[i]);
		}
		for (; i < b.size; i++) {
			carry = ADDCARRY(carry, b.data[i], 0, &out[i]);
		}
		if (assign_carry && carry) out[b.size] = carry;
		return b.size + carry;
	}
}

// out = a - b;
// a must be >= b
static size_t usub(Slice a, Slice b, DataBlock* out) {
	assert(cmp(a, b) >= 0);
	assert(a.size == size_without_zeros(a.data, a.size));
	assert(b.size == size_without_zeros(b.data, b.size));
	unsigned char borrow = 0;
	size_t i = 0;
	for (; i < b.size; i++) {
		borrow = SUBBORROW(borrow, a.data[i], b.data[i], &out[i]);
	}
	for (; i < a.size; i++) {
		borrow = SUBBORROW(borrow, a.data[i], 0, &out[i]);
	}
	return size_without_zeros(out, a.size);
}

// usub_nlz - unsigned subtraction, no leading zeros;
// out = a - b;
// does not assign leading zeros;
// a must be >= b
static size_t usub_nlz(Slice a, Slice b, DataBlock* out) {
	assert(cmp(a, b) >= 0);
	assert(a.size == size_without_zeros(a.data, a.size));
	assert(b.size == size_without_zeros(b.data, b.size));
	unsigned char borrow = 0;
	size_t i = 0;
	size_t zeros = 0;
	for (; i < b.size; i++) {
		borrow = subborrow_nlz(borrow, a.data[i], b.data[i], &out[i], &zeros);
	}
	for (; i < a.size; i++) {
		borrow = subborrow_nlz(borrow, a.data[i], 0, &out[i], &zeros);
	}
	return i - zeros;
}

// out = a - b;
// a and b are unsigned but output can be negative;
// nlz (no leading zeros): if true uses usub_nlz, if false uses usub
static size_t sub(Slice a, Slice b, DataBlock* out, bool* out_sign, bool nlz) {
	size_t diff_point;
	int cmp_res = cmpd(a, b, &diff_point);
	if (cmp_res == 0) {
		if (out_sign) *out_sign = 0;
		return 0;
	}
	split(&a, diff_point);
	split(&b, diff_point);
	if (a.size == 0) {
		*out_sign = 1;
		memmove(out, b.data, b.size * sizeof(DataBlock));
		return b.size;
	}
	if (b.size == 0) {
		*out_sign = 0;
		memmove(out, a.data, a.size * sizeof(DataBlock));
		return a.size;
	}
	if (cmp_res > 0) {
		if (out_sign) *out_sign = 0;
		size_t size = nlz ? usub_nlz(a, b, out) : usub(a, b, out);
		assert(size == size_without_zeros(out, size));
		return size;
	} else {
		if (out_sign) *out_sign = 1;
		size_t size = nlz ? usub_nlz(b, a, out) : usub(b, a, out);
		assert(size == size_without_zeros(out, size));
		return size;
	}
}

// out = (a << lshift) - b;
// Requirements: (a << lshift) >= b;
// NOTE: lshift by number of blocks not bits
static size_t usub_lshifted(Slice a, size_t lshift, Slice b, DataBlock* out) {
	assert(cmp_lshifted(a, lshift, b) >= 0);
	unsigned char borrow = 0;
	size_t i = 0;
	if (b.size <= lshift) {
		for (; i < b.size; i++) {
			borrow = SUBBORROW(borrow, 0, b.data[i], &out[i]);
		}
		for (; i < lshift; i++) {
			borrow = SUBBORROW(borrow, 0, 0, &out[i]);
		}
	} else {
		for (; i < lshift; i++) {
			borrow = SUBBORROW(borrow, 0, b.data[i], &out[i]);
		}
		for (; i < b.size; i++) {
			borrow = SUBBORROW(borrow, a.data[i - lshift], b.data[i], &out[i]);
		}
	}
	for (; i < a.size + lshift; i++) {
		borrow = SUBBORROW(borrow, a.data[i - lshift], 0, &out[i]);
	}
	return size_without_zeros(out, i);
}

// out = (a << lshift) - b;
// assigns no leading zeros;
// Requirements: (a << lshift) >= b;
// NOTE: lshift by number of blocks not bits
static size_t usub_lshifted_nlz(Slice a, size_t lshift, Slice b, DataBlock* out) {
	assert(cmp_lshifted(a, lshift, b) >= 0);
	unsigned char borrow = 0;
	size_t i = 0;
	size_t zeros = 0;
	if (b.size <= lshift) {
		for (; i < b.size; i++) {
			borrow = subborrow_nlz(borrow, 0, b.data[i], &out[i], &zeros);
		}
		for (; i < lshift; i++) {
			borrow = subborrow_nlz(borrow, 0, 0, &out[i], &zeros);
		}
	} else {
		for (; i < lshift; i++) {
			borrow = subborrow_nlz(borrow, 0, b.data[i], &out[i], &zeros);
		}
		for (; i < b.size; i++) {
			borrow = subborrow_nlz(borrow, a.data[i - lshift], b.data[i], &out[i], &zeros);
		}
	}
	for (; i < a.size + lshift; i++) {
		borrow = subborrow_nlz(borrow, a.data[i - lshift], 0, &out[i], &zeros);
	}
	return i - zeros;
}

// out = a - (b << lshift);
// a must be >= (b << lshift);
// NOTE: lshift by number of blocks not bits
static size_t usub_b_lshifted(Slice a, Slice b, size_t lshift, DataBlock* out) {
	assert(cmp_lshifted(b, lshift, a) <= 0);
	memmove(out, a.data, lshift * sizeof(DataBlock));
	unsigned char borrow = 0;
	size_t i = lshift;
	for (; i < b.size + lshift; i++) {
		borrow = SUBBORROW(borrow, a.data[i], b.data[i - lshift], &out[i]);
	}
	for (; i < a.size; i++) {
		borrow = SUBBORROW(borrow, a.data[i], 0, &out[i]);
	}
	return size_without_zeros(out, i);
}

// out = a - (b << lshift);
// assigns no leading zeros;
// a must be >= (b << lshift);
// NOTE: lshift by number of blocks not bits
static size_t usub_b_lshifted_nlz(Slice a, Slice b, size_t lshift, DataBlock* out) {
	assert(cmp_lshifted(b, lshift, a) <= 0);
	memmove(out, a.data, lshift * sizeof(DataBlock));
	unsigned char borrow = 0;
	size_t i = lshift;
	size_t zeros = 0;
	for (; i < b.size + lshift; i++) {
		borrow = subborrow_nlz(borrow, a.data[i], b.data[i - lshift], &out[i], &zeros);
	}
	for (; i < a.size; i++) {
		borrow = subborrow_nlz(borrow, a.data[i], 0, &out[i], &zeros);
	}
	return i - zeros;
}

static size_t sub_lshifted(Slice a, size_t lshift, Slice b, DataBlock* out, bool* out_sign, bool nlz) {
	size_t diff_point;
	int cmp_res = cmpd_lshifted(a, lshift, b, &diff_point);
	if (cmp_res == 0) {
		if (out_sign) *out_sign = 0;
		return 0;
	}
	split(&b, diff_point);
	if (diff_point > lshift) {
		split(&a, diff_point - lshift);
	} else {
		a.size = 0;
	}
	if (a.size == 0) {
		*out_sign = 1;
		memmove(out, b.data, b.size * sizeof(DataBlock));
		return b.size;
	}
	if (b.size == 0) {
		*out_sign = 0;
		memset(out, 0, lshift * sizeof(DataBlock));
		memmove(out + lshift, a.data, a.size * sizeof(DataBlock));
		return a.size + lshift;
	}
	if (cmp_res > 0) {
		if (out_sign) *out_sign = 0;
		size_t size = nlz
			? usub_lshifted_nlz(a, lshift, b, out)
			: usub_lshifted(a, lshift, b, out);
		assert(size == size_without_zeros(out, size));
		return size;
	} else {
		if (out_sign) *out_sign = 1;
		size_t size = nlz
			? usub_b_lshifted_nlz(b, a, lshift, out)
		    : usub_b_lshifted(b, a, lshift, out);
		assert(size == size_without_zeros(out, size));
		return size;
	}
}

static size_t add_signed(Slice a, bool a_sign, Slice b, bool b_sign,
		DataBlock* out, bool* out_sign, bool assign_carry, bool nlz) {
	if (a_sign == b_sign) {
		*out_sign = a_sign;
		return add(a, b, out, assign_carry);
	}
	size_t size = sub(a, b, out, out_sign, nlz);
	*out_sign ^= a_sign;
	return size;
}

// (a << lshift) + b, signed
// NOTE: lshift by number of blocks not bits
static size_t add_lshifted_signed(Slice a, size_t lshift, bool a_sign,
		Slice b, bool b_sign, DataBlock* out, bool* out_sign, bool assign_carry, bool nlz) {
	if (a_sign == b_sign) {
		*out_sign = a_sign;
		return add_lshifted(a, lshift, b, out, assign_carry);
	}
	size_t size = sub_lshifted(a, lshift, b, out, out_sign, nlz);
	*out_sign ^= a_sign;
	return size;
}

static size_t mul_basic(Slice a, Slice b, DataBlock* out) {
	if (a.size == 0 || b.size == 0) return 0;
	const size_t max_size = a.size + b.size;
	memset(out, 0, max_size * sizeof(DataBlock));
	for (size_t i = 0; i < a.size; i++) {
		DataBlock d1 = a.data[i];
		for (size_t j = 0; j < b.size; j++) {
			DataBlock d2 = b.data[j];
			DataBlock high;
			DataBlock low = block_mul(d1, d2, &high);
			unsigned char carry;
			carry = ADDCARRY(0, out[i + j], low, &out[i + j]);
			carry = ADDCARRY(carry, out[i + j + 1], high, &out[i + j + 1]);
			size_t idx = i + j + 2;
			while (carry) {
				carry = ADDCARRY(carry, out[idx], 0, &out[idx]);
				idx++;
			}
		}
	}
	return size_without_zeros(out, max_size);
}


static size_t karatsuba(Slice b, Slice d, size_t n, DataBlock* out, DataBlock* buffer);

static size_t mul(Slice a, Slice b, DataBlock* out, DataBlock* buffer) {
	assert(size_without_zeros(a.data, a.size) == a.size);
	assert(size_without_zeros(b.data, b.size) == b.size);
	if (a.size == 0 || b.size == 0) return 0;
	const int n = a.size > b.size ? a.size : b.size;
	if (n >= BIGINT_KARATSUBA_THRESHOLD) {
/*
		for (size_t i = a.size; i > 0; i--) {
			printf(FORMAT_HEX " ", b.data[i - 1]);
		} printf ("* ");
		for (size_t i = b.size; i > 0; i--) {
			printf(FORMAT_HEX " ", b.data[i - 1]);
		} printf ("= ");
*/
		size_t size = karatsuba(a, b, n, out, buffer);
		assert(size_without_zeros(out, size) == size);
/*
		for (size_t i = size; i > 0; i--) {
			printf(FORMAT_HEX " ", out[i - 1]);
		} printf ("\n");
*/
		return size;
	}
	return mul_basic(a, b, out);
}

// returns buffer size required for karatsuba in bytes
static inline size_t karatsuba_buffer_size(size_t n) {
	size_t sum = 0;
	constexpr size_t threshold = 2 > BIGINT_KARATSUBA_THRESHOLD
		? 2 : BIGINT_KARATSUBA_THRESHOLD;
	while (n >= threshold) {
		if (n & 1) n += 1;
		sum += n;
		n >>= 1;
	}
	return sum * sizeof(DataBlock);
}

// Karatsuba multiplication
// @params:
// b - first number (or slice of a number);
// d - second number (or slice);
// n - max(b.size, d.size);
// out - filled with b * d, must be allocated outside
//       (at least (b.size + d.size) * sizeof(DataBlock) bytes);
// buffer - memory for storing intermediate calculations
//          (karatsuba_buffer_size() * sizeof(DataBlock) bytes)
// @return size of out
//
// Note: a slightly modified formula -
// (a-b)(d-c)+ac+bd is used for the middle part because
// (a+b) and (c+d) can overflow
static size_t karatsuba(Slice b, Slice d, size_t n, DataBlock* out, DataBlock* buffer) {
	if(b.size == 0 || d.size == 0) return 0;
	assert(b.size == size_without_zeros(b.data, b.size));
	assert(d.size == size_without_zeros(d.data, d.size));
	assert(n);

	if (n == 1) {
		DataBlock high = 0;
		out[0] = block_mul(b.data[0], d.data[0], &high);
		assert(high || out[0]);
		if (high) {
			out[1] = high;
			return 2;
		} else {
			return 1;
		}
	}

	const size_t m = (n + 1) / 2;
	assert(m);
	assert(m <= n);
	const Slice a = split(&b, m);
	const Slice c = split(&d, m);

	const size_t ab_cap = m;
	const size_t dc_cap = m;
	const size_t abcd_cap = 2 * m;
	const size_t ac_cap = 2 * (n - m);
	const size_t bd_cap = 2 * m;

	DataBlock* const abcd_ptr = out;
	DataBlock* const ab_ptr = buffer;
	DataBlock* const dc_ptr = buffer + m;
	DataBlock* const ac_ptr = buffer;
	DataBlock* const bd_ptr = buffer;
	DataBlock* const buf_ptr = buffer + 2 * m;

	bool ab_sign;
	bool dc_sign;
	bool abcd_sign;
	
	const Slice a_minus_b = {
		.data = ab_ptr,
		.size = sub(a, b, ab_ptr, &ab_sign, false)
	};
	assert(a_minus_b.size <= ab_cap);
	assert(a_minus_b.size == size_without_zeros(a_minus_b.data, a_minus_b.size));

	const Slice d_minus_c = {
		.data = dc_ptr,
		.size = sub(d, c, dc_ptr, &dc_sign, false)
	};
	assert(d_minus_c.size <= dc_cap);
	assert(d_minus_c.size == size_without_zeros(d_minus_c.data, d_minus_c.size));

	abcd_sign = ab_sign ^ dc_sign;

	// abcd = (a-b)(d-c)
	Slice abcd = {
		.data = abcd_ptr,
		.size = mul(a_minus_b, d_minus_c, abcd_ptr, buf_ptr)
	};
	assert(abcd.size <= abcd_cap);
	assert(abcd.size == size_without_zeros(abcd.data, abcd.size));

	const Slice a_times_c = { 
		.data = ac_ptr,
		.size = mul(a, c, ac_ptr, buf_ptr)
	};
	assert(a_times_c.size <= ac_cap);
	assert(a_times_c.size == size_without_zeros(a_times_c.data, a_times_c.size));

	// abcd = abcd + ac
	abcd.size = add_signed(abcd, abcd_sign, a_times_c, 0, abcd_ptr, &abcd_sign, true, true);
	assert(abcd.size == size_without_zeros(abcd.data, abcd.size));
	// abcd = (a_times_c << m) + abcd
	abcd.size = add_lshifted_signed(a_times_c, m, 0, abcd, abcd_sign, abcd_ptr, &abcd_sign, true, true);
	assert(abcd.size == size_without_zeros(abcd.data, abcd.size));

	const Slice b_times_d = {
		.data = bd_ptr,
		.size = mul(b, d, bd_ptr, buf_ptr)
	};
	assert(b_times_d.size <= bd_cap);
	assert(b_times_d.size == size_without_zeros(b_times_d.data, b_times_d.size));

	// abcd = abcd + bd
	abcd.size = add_signed(abcd, abcd_sign, b_times_d, 0, abcd_ptr, &abcd_sign, true, true);
	assert(abcd.size == size_without_zeros(abcd.data, abcd.size));
	// move abcd from out to out + m
	memmove(out + m, out, abcd.size * sizeof(DataBlock));
	abcd.data += m;
	// (abcd << m) + bd
	return add_lshifted_signed(abcd, m, abcd_sign, b_times_d, 0, out, &abcd_sign, true, true);
}

// |out| = |a| + |b|
BigInt bigint_uadd(ConstBigInt a, ConstBigInt b, BigInt* out_ptr) {
	assert(a);
	assert(b);
	assert(out_ptr);

	// add_a() needs a.size >= b.size, so swap a and b if otherwise
	if (a->size < b->size) {
		ConstBigInt tmp = a;
		a = b;
		b = tmp;
	}

	BigInt out = bigint_realloc_if_small(out_ptr, a->size, false);
	if (!out) return NULL;

	Slice A = { .data = a->data, .size = a->size };
	Slice B = { .data = b->data, .size = b->size };

	out->size = add_a(A, B, out->data, false);

	// carry happened
	if (out->size > A.size) {
		out = bigint_realloc_if_small(out_ptr, A.size + 1, true);
		if (!out) return NULL;
		out->data[A.size] = 1;
	}

	return out;
}

// |out| = |a| - |b|,
// assumes |a| >= |b|
BigInt bigint_usub(ConstBigInt a, ConstBigInt b, BigInt* out_ptr) {
	assert(a);
	assert(b);
	assert(out_ptr);

	Slice A = { .data = a->data, .size = a->size };
	Slice B = { .data = b->data, .size = b->size };
	int cmp_res = cmp(A, B);
	if (cmp_res == 0) {
		return bigint_set_zero(out_ptr);
	}
	if (cmp_res < 0) {
		bigint_errno = BIGINT_ERR_SUB_NEG;
		return NULL;
	}
	BigInt out;
	if (!(out = bigint_realloc_if_small(out_ptr, a->size, false))) return NULL;
	bool sign;
	out->size = sub(A, B, out->data, &sign, false);
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

	if (a->size == 0 || b->size == 0) {
		return bigint_set_zero(out_ptr);
	}

	const size_t n = (a->size > b->size) ? a->size : b->size;
	const size_t cap = a->size + b->size;

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

	out = bigint_realloc_if_small(&out, cap, false);
	if (!out) goto error;

	Slice A = { .data = a->data, .size = a->size };
	Slice B = { .data = b->data, .size = b->size };
	size_t bufsize = karatsuba_buffer_size(n);
	if (bufsize > ctx.bufsize) {
		ctx.buf = realloc(ctx.buf, bufsize);
	}
	out->size = karatsuba(A, B, n, out->data, ctx.buf);
	assert(out->size <= cap);
	assert(out->size == size_without_zeros(out->data, out->size));
	/*
	BigInt check = bigint_alloc(cap);
	check->size = mul_basic(A, B, check->data);
	assert(bigint_cmp(out, check) == 0 || ({
				bigint_printf("%px, %px\n", out, check);
				0;
			}));
	bigint_free(check);
	*/

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
		if (!bigint_lshift_blocks(b, size_diff, &c)) goto error;

		if (bigint_ucmp(rem, c) < 0) {
			if (!bigint_rezalloc(&out, size_diff)) goto error;
			out->size = size_diff;
			bigint_rshift_blocks(c, 1, &c);
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

BigInt bigint_lshift_blocks(ConstBigInt z, size_t shift, BigInt* out_ptr) {
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

BigInt bigint_rshift_blocks(ConstBigInt z, size_t shift, BigInt* out_ptr) {
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

	const size_t shift_blocks = shift / BLOCK_WIDTH;
	const size_t shift_bits = shift % BLOCK_WIDTH;
	if (shift_bits == 0) return bigint_lshift_blocks(z, shift_blocks, out_ptr);

	const DataBlock LOW_BITS = (1ULL << (BLOCK_WIDTH - shift_bits)) - 1;
	const DataBlock HIGH_BITS = ~0ULL - LOW_BITS;

	size_t out_size = z->size + shift_blocks;

	size_t i = z->size - 1;
	DataBlock hi, lo;

	hi = z->data[i] & HIGH_BITS;
	lo = z->data[i] & LOW_BITS;
	BigInt out;
	if (hi) {
		out = bigint_resize(out_ptr, ++out_size);
		if (!out) return NULL;
		out->data[i + shift_blocks + 1] = hi >> (BLOCK_WIDTH - shift_bits);
	} else {
		out = bigint_resize(out_ptr, out_size);
		if (!out) return NULL;
	}
	out->data[i + shift_blocks] = lo << shift_bits;

	if (i--) while (1) {
		hi = z->data[i] & HIGH_BITS;
		lo = z->data[i] & LOW_BITS;
		out->data[i + shift_blocks] = lo << shift_bits;
		out->data[i + shift_blocks + 1] |= hi >> (BLOCK_WIDTH - shift_bits);
		if (i == 0) break;
		i--;
	}

	memset(out->data, 0, shift_blocks * sizeof(DataBlock));
	return out;
}

BigInt bigint_rshift(ConstBigInt z, size_t shift, BigInt* out_ptr) {
	assert(z);
	assert(out_ptr);

	size_t shift_blocks = shift / BLOCK_WIDTH;
	size_t shift_bits = shift % BLOCK_WIDTH;
	if (shift_bits == 0) return bigint_rshift_blocks(z, shift_blocks, out_ptr);

	size_t z_size = z->size;
	if (shift_blocks >= z_size) {
		return bigint_set_zero(out_ptr);
	}

	size_t out_size = z_size - shift_blocks;
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
		hi = z->data[i + shift_blocks] & HIGH_BITS;
		lo = z->data[i + shift_blocks + 1] & LOW_BITS;
		out->data[i] = (hi >> shift_bits) | lo << (BLOCK_WIDTH - shift_bits);
	}

	hi = z->data[i + shift_blocks] & HIGH_BITS;
	if (i + shift_blocks + 1 < z_size)
		lo = z->data[i + shift_blocks + 1] & LOW_BITS;
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
					bifs.leading_zeros = true;
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
	if (bifs.leading_zeros) offset += sprintf(str + offset, f_0hex, z->data[i]);
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
