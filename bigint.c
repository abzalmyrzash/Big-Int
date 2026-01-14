#include "bigint_impl.h"
#include "elog.h"
#include <assert.h>
#include <errno.h>
#include <stdlib.h>
#include <stdarg.h>
#include <time.h>

int bigint_errno = 0;

static unsigned long long cnt_alloc = 0;
static unsigned long long cnt_realloc = 0;
static unsigned long long cnt_free = 0;

#define AS_SLICE(x) (Slice) _Generic((x),                         \
	BigInt:      (Slice){ .data = (x)->data, .size = (x)->size }, \
	ConstBigInt: (Slice){ .data = (x)->data, .size = (x)->size })

void bigint_init(void) {
	arena = arena_create(GiB(1), sizeof(mem_arena));
	int i = arena->pos;

	INITIAL_NEWTON_PRECISION = log2(17);
	NEWTON_STEPS = 0;

	_48_OVER_17.data = malloc(sizeof(Block));
	_32_OVER_17.data = malloc(sizeof(Block));
	_48_OVER_17_REM.data = malloc(sizeof(Block));
	_32_OVER_17_REM.data = malloc(sizeof(Block));

	_48_OVER_17_POINT = INITIAL_NEWTON_PRECISION;
	_48_OVER_17_WIDTH = _48_OVER_17_INT_WIDTH + _48_OVER_17_POINT;
	MUT_DATA(_48_OVER_17)[0] = (48ULL << _48_OVER_17_POINT) / 17;
	_48_OVER_17.size = 1;
	MUT_DATA(_48_OVER_17_REM)[0] = (48ULL << _48_OVER_17_POINT) % 17;
	_48_OVER_17_REM.size = 1;

	_32_OVER_17_POINT = INITIAL_NEWTON_PRECISION;
	_32_OVER_17_WIDTH = _32_OVER_17_INT_WIDTH + _32_OVER_17_POINT;
	MUT_DATA(_32_OVER_17)[0] = (32ULL << _32_OVER_17_POINT) / 17;
	_32_OVER_17.size = 1;
	MUT_DATA(_32_OVER_17_REM)[0] = (32ULL << _32_OVER_17_POINT) % 17;
	_32_OVER_17_REM.size = 1;

	/*
	print_slice_bin_point(_48_OVER_17, _48_OVER_17_POINT);
	printf(":");
	print_slice_bin(_48_OVER_17_REM);
	printf("\n");
	print_slice_bin_point(_32_OVER_17, _32_OVER_17_POINT);
	printf(":");
	print_slice_bin(_32_OVER_17_REM);
	printf("\n");
	*/

}

void bigint_finish(void) {
	if (arena->pos != ARENA_BASE_POS) {
		WLOG("WARNING: Arena hasn't been fully cleared (pos = %zu)\n", arena->pos);
	}
	arena_destroy(arena);
}

BigInt bigint_alloc(size_t cap) {
	if (cap > MAX_CAP) {
		ELOG_STR("cap exceeds MAX_CAP");
		ELOG("Parameters: new_size = %zu\n", cap);
		return NULL;
	}
	BigInt z = malloc(DATA_OFFSET + cap * sizeof(Block));
	if (z == NULL) {
		bigint_errno = BIGINT_ERR_ALLOC_FAIL;
		ELOG_CODE(errno);
		ELOG("Parameters: cap = %zu\n", cap);
		return NULL;
	}
	cnt_alloc++;
	z->cap = cap;
	z->size = 0;
	z->sign = 0;
	return z;
}

BigInt bigint_zalloc(size_t cap) {
	if (cap > MAX_CAP) {
		ELOG_STR("cap exceeds MAX_CAP");
		ELOG("Parameters: new_size = %zu\n", cap);
		return NULL;
	}
	BigInt z = calloc(1, DATA_OFFSET + cap * sizeof(Block));
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
		z = realloc(z, DATA_OFFSET + cap * sizeof(Block));
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

BigInt bigint_reserve(BigInt* z_ptr, size_t cap, bool keep_data) {
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
		memset(z->data, 0, cap * sizeof(Block));
		z->size = 0;
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
	z = bigint_reserve(z_ptr, new_size, KEEP_DATA);
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

void bigint_structinfo(void) {
	printf("sizeof(bigint): %zu\n", sizeof(struct bigint));
	printf("offsetof(bigint, data): %zu\n", offsetof(struct bigint, data));
	printf("sizeof(Block): %zu\n", sizeof(Block));
}

void bigint_memstat(void) {
	printf("Alloc: %llu\n", cnt_alloc);
	printf("Realloc: %llu\n", cnt_realloc);
	printf("Free: %llu\n", cnt_free);
}

size_t bigint_size(ConstBigInt z) {
	return z->size;
}

size_t bigint_width(ConstBigInt z) {
	if (z->size == 0) return 0;
	return z->size * BLOCK_WIDTH - CLZ(z->data[z->size - 1]);
}

size_t bigint_cap(ConstBigInt z) {
	return z->cap;
}

const Block* bigint_data(ConstBigInt z) {
	return z->data;
}

bool bigint_sign(ConstBigInt z) {
	return z->sign;
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
	z->sign = (v < 0);
	return z;
}

BigInt bigint_set_usmall(BigInt* z_ptr, USmallInt v) {
	if (v == 0) return bigint_set_zero(z_ptr);
	BigInt z = bigint_resize(z_ptr, 1);
	if (!z) return NULL;
	z->data[0] = v;
	z->sign = 0;
	return z;
}

BigInt bigint_create_small(SmallInt v) {
	if (v == 0) return bigint_create_zero();
	BigInt z = bigint_alloc(1);
	z->data[0] = ABS(v);
	z->size = 1;
	z->sign = (v < 0);
    return z;
}

BigInt bigint_create_usmall(USmallInt v) {
	if (v == 0) return bigint_create_zero();
	BigInt z = bigint_alloc(1);
	z->data[0] = v;
	z->size = 1;
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
	memcpy(dst->data, src->data, src->size * sizeof(Block));
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

	BigInt out = bigint_reserve(out_ptr, a->size, CLEAR_DATA);
	if (!out) return NULL;

	Slice A = { .data = a->data, .size = a->size };
	Slice B = { .data = b->data, .size = b->size };

	out->size = add_a(A, B, out->data, false);

	// carry happened
	if (out->size > A.size) {
		out = bigint_reserve(out_ptr, A.size + 1, KEEP_DATA);
		if (!out) return NULL;
		out->data[A.size] = 1;
	}

	return out;
}

// out = |a| - |b|
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
	BigInt out = bigint_reserve(out_ptr, a->size, CLEAR_DATA);
	if (!out) return NULL;
	bool sign;
	out->size = sub(A, B, out->data, &sign);
	out->sign = sign;
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
			bigint_set_zero(out_ptr);
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
			bigint_set_zero(out_ptr);
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

	Slice A = AS_SLICE(a);
	Slice B = AS_SLICE(b);

	u64 restore_pos = arena->pos;
	if (out == a) {
		Block* tmp = PUSH_ARRAY(arena, Block, A.size);
		A.data = memcpy(tmp, A.data, A.size * sizeof(Block));
	}
	if (out == b) {
		Block* tmp = PUSH_ARRAY(arena, Block, B.size);
		B.data = memcpy(tmp, B.data, B.size * sizeof(Block));
	}

	out = bigint_reserve(&out, cap, CLEAR_DATA);
	if (!out) goto error;

	out->size = mul(A, B, out->data);
	assert(out->size <= cap);
	assert(out->size == size_without_zeros(out->data, out->size));
	assert(({
		Block* check = malloc(cap * sizeof(Block));
		size_t check_size = long_mul(A, B, check);
		int cmp_res = cmp((Slice){out->data, out->size}, (Slice){check, check_size});
		free(check);
		cmp_res == 0;
	}));

	arena_pop_to(arena, restore_pos);

	return *out_ptr = out;
error:
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

	u64 restore_pos = arena->pos;

	BigInt out = *out_ptr;
	BigInt rem = rem_ptr ? *rem_ptr : NULL;

	Slice A = { .data = a->data, .size = a->size };
	Slice B = { .data = b->data, .size = b->size };

	Block* rem_data;
	size_t rem_size;

	if (out == a) {
		Block* tmp = PUSH_ARRAY(arena, Block, A.size);
		A.data = memcpy(tmp, A.data, A.size * sizeof(Block));
	} else if (out == b) {
		Block* tmp = PUSH_ARRAY(arena, Block, B.size);
		B.data = memcpy(tmp, B.data, B.size * sizeof(Block));
	}

	out = bigint_reserve(&out, A.size - B.size + 1, CLEAR_DATA);
	if (!out) goto error;

	if (rem_ptr) {
		rem = bigint_reserve(&rem, A.size, CLEAR_DATA);
		if (!rem) goto error;
		rem_data = rem->data;
		if (rem == b) {
			Block* tmp = PUSH_ARRAY(arena, Block, B.size);
			B.data = memcpy(tmp, B.data, B.size * sizeof(Block));
		}
	} else {
		rem_data = PUSH_ARRAY(arena, Block, A.size);
	}

	out->size = long_div(A, B, out->data, rem_data, &rem_size);
	if (rem) rem->size = rem_size;

ret:
	arena_pop_to(arena, restore_pos);
	if (rem_ptr) *rem_ptr = rem;
	return *out_ptr = out;

error:
	arena_pop_to(arena, restore_pos);
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
	memmove((*out_ptr)->data + shift, z->data, z->size * sizeof(Block));
	memset((*out_ptr)->data, 0, shift * sizeof(Block));
	return *out_ptr;
}

BigInt bigint_rshift_blocks(ConstBigInt z, size_t shift, BigInt* out_ptr) {
	assert(z);
	assert(out_ptr);
	if (z->size <= shift) return bigint_set_zero(out_ptr);
	if (shift == 0) return bigint_copy(out_ptr, z);

	size_t out_size = z->size - shift;
	if (!bigint_resize(out_ptr, out_size)) return NULL;
	memmove((*out_ptr)->data, z->data + shift, out_size * sizeof(Block));
	return *out_ptr;
}

BigInt bigint_lshift(ConstBigInt z, size_t shift, BigInt* out_ptr) {
	assert(z);
	assert(out_ptr);
	if (z->size == 0) return bigint_set_zero(out_ptr);

	const size_t shift_blocks = shift / BLOCK_WIDTH;
	const size_t shift_bits = shift % BLOCK_WIDTH;
	if (shift_bits == 0) return bigint_lshift_blocks(z, shift_blocks, out_ptr);

	const Block LOW_BITS = (1ULL << (BLOCK_WIDTH - shift_bits)) - 1;
	const Block HIGH_BITS = ~0ULL - LOW_BITS;

	size_t out_size = z->size + shift_blocks;

	size_t i = z->size - 1;
	Block hi, lo;

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

	memset(out->data, 0, shift_blocks * sizeof(Block));
	return out;
}

BigInt bigint_rshift(ConstBigInt z, size_t shift, BigInt* out_ptr) {
	assert(z);
	assert(out_ptr);

	const size_t shift_blocks = shift / BLOCK_WIDTH;
	const int shift_bits = shift % BLOCK_WIDTH;
	if (shift_bits == 0) return bigint_rshift_blocks(z, shift_blocks, out_ptr);

	size_t z_size = z->size;
	if (shift_blocks >= z_size) {
		return bigint_set_zero(out_ptr);
	}

	size_t out_size = z_size - shift_blocks;
	if (shift_bits >= BLOCK_WIDTH - CLZ((Block)z->data[z_size - 1])) {
		if (out_size <= 1) return bigint_set_zero(out_ptr);
		out_size--;
	}
	BigInt out = bigint_resize(out_ptr, out_size);
	if (!out) return NULL;

	const Block LOW_BITS = (1ULL << shift_bits) - 1;
	const Block HIGH_BITS = ~0ULL - LOW_BITS;

	size_t i;
	Block hi, lo;
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

BigInt bigint_pow(ConstBigInt a, size_t exp, BigInt* out_ptr) {
	assert(a);
	assert(exp);
	assert(out_ptr);
	
	BigInt out = bigint_reserve(out_ptr, to_blocks(exp * bigint_width(a)), CLEAR_DATA);
	if (!out) return NULL;

	out->size = power(AS_SLICE(a), exp, out->data);
	return out;
}

BigInt bigint_powmod(ConstBigInt a, ConstBigInt exp, ConstBigInt mod, BigInt* out_ptr) {
	assert(a);
	assert(exp);
	assert(mod);
	assert(out_ptr);
	
	BigInt out = bigint_reserve(out_ptr, mod->size, CLEAR_DATA);
	if (!out) return NULL;

	out->size = powmod(AS_SLICE(a), AS_SLICE(exp), AS_SLICE(mod), out->data);
	return out;
}

BigInt bigint_powmod_recpr(ConstBigInt a, ConstBigInt exp, ConstBigInt mod, ConstBigInt recpr, size_t precision, BigInt* out_ptr) {
	assert(a);
	assert(exp);
	assert(mod);
	assert(out_ptr);
	
	BigInt out = bigint_reserve(out_ptr, mod->size, CLEAR_DATA);
	if (!out) return NULL;

	out->size = powmod_recpr(AS_SLICE(a), AS_SLICE(exp), AS_SLICE(mod), AS_SLICE(recpr), precision, out->data);
	return out;
}

int bigint_fwrite(FILE* file, ConstBigInt z, bool is_signed) {
	assert(z);
	uint64_t size = z->size;
	if (fwrite(&size, sizeof(size), 1, file) != 1) {
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

	bool sign;
	if (is_signed && fread(&sign, sizeof(sign), 1, file) != 1) return errno;
	z->sign = sign;

	if (fread(&z->data, sizeof(*z->data), size, file) != size) return errno;
	return 0;
}

FormatSpec bigint_initBIFS(void) {
	return (FormatSpec) { 0 };
}

static int _bigint_fprintf(FILE* file, const char* const format, va_list args) {
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
					res = putc('%', file);
					if (res < 0) {
						ELOG_CODE_STR(errno, "putc failed");
						return res;
					}
					processing_bifs = false;
					bifs = bigint_initBIFS();
					break;

				default:
					ELOG_STR("Invalid format specifier '%c' in \"%s\".", c, format);
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
				ELOG_CODE_STR(res, "putc failed.");
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
		(z->size > 0) ? bigint_decimal_width(z) : 1;
	const size_t max_spaces = (max_digits / DEC_DIGITS_PER_SPACE) * bifs.add_spaces;
	const size_t max_len = max_digits + max_spaces + add_sign;
	char* const str = malloc(max_len + 1);

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

	char* digits = str + offset;

	size_t len = decimal_split(AS_SLICE(z), digits);
	digits[len] = '\0';

	// 1234567890123
	// 1    234567890123
	// 1 234   567890123
	// 1 234 567  890123
	// 1 234 567 890 123

	if (bifs.add_spaces) {
		size_t spaces = len / DDPS;
		size_t space1 = ((len - 1) % DDPS) + 1;
		digits += space1;
		len -= space1;
		memmove(digits + spaces, digits, sizeof(char) * len);
		digits[0] = ' ';
		for (size_t i = 0; i < spaces - 1; i++) {
			memmove(digits + 1 + i * (DDPS + 1), digits + spaces + i * DDPS, sizeof(char) * DDPS);
			digits[(i + 1) * (DDPS + 1)] = ' ';
		}
		digits[spaces * (DDPS + 1)] = '\0';
	}
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
	if (str_len == 0) return bigint_set_zero(out_ptr);
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
				bit_set(out->data, shift);
			if (value & 2)
				bit_set(out->data, shift + 1);
			if (value & 4)
				bit_set(out->data, shift + 2);
			if (value & 8)
				bit_set(out->data, shift + 3);
		}
		shift += 4;
	}

	return out;
}

BigInt bigint_scan_dec(const char* str, size_t str_len, BigInt* out_ptr) {
	assert(str);
	assert(out_ptr);
	if (str_len == 0) return bigint_set_zero(out_ptr);
	size_t cap = str_len / (MAX_DEC_DIGITS_PER_BLOCK - 1) + 1;
	assert(cap > 0);
	BigInt out = bigint_rezalloc(out_ptr, cap);
	if (!out) return NULL;

	const size_t pow10_cap = cap + 1;
	const size_t digit_cap = cap + 1;
	const size_t tmp_cap = cap;

	u64 restore_pos = arena->pos;
	Block* const pow10_data = PUSH_ARRAY(arena, Block, pow10_cap);
	Block* const digit_data = PUSH_ARRAY(arena, Block, digit_cap);
	Block* const tmp_data = PUSH_ARRAY(arena, Block, tmp_cap);
	Slice pow10 = { .data = pow10_data, .size = 1 };
	pow10_data[0] = 1;
	Slice digit = { .data = digit_data };
	Slice tmp = { .data = tmp_data };

	size_t i = str_len - 1;
	while (1) {
		const char ch = str[i];
		if (ch < '0' || ch > '9') {
			i--;
			continue;
		}
		tmp.size = set_usmall(tmp_data, ch - '0');

		digit.size = mul(tmp, pow10, digit_data);
		out->size = add((Slice){out->data, out->size}, digit, out->data, true);
		if (i == 0) break;
		i--;
		tmp.size = move(pow10, tmp_data);
		pow10.size = mul(tmp, TEN, pow10_data);
	}

	arena_pop_to(arena, restore_pos);
	return out;
}

BigInt bigint_urecpr(ConstBigInt d, size_t precision, BigInt* out_ptr) {
	assert(d);
	assert(out_ptr);
	if (d->size == 0) {
		bigint_errno = BIGINT_ERR_DIV_BY_ZERO;
		ELOG_STR("DIVISION BY ZERO");
		return NULL;
	}
	size_t size = newton_reciprocal_cap(precision);
	BigInt out = bigint_reserve(out_ptr, size, CLEAR_DATA);
	if (!out) return NULL;
	out->size = newton_reciprocal(AS_SLICE(d), precision, out->data);
	return out;
}

BigInt bigint_recpr(ConstBigInt d, size_t precision, BigInt* out_ptr) {
	assert(d);
	assert(out_ptr);
	BigInt out = bigint_urecpr(d, precision, out_ptr);
	if (!out) return NULL;
	out->sign = d->sign;
	return out;
}

BigInt bigint_udiv_recpr(ConstBigInt num, ConstBigInt denum, ConstBigInt recpr, size_t precision, BigInt* quo, BigInt* rem) {
	assert(num);
	assert(denum);
	assert(recpr);
	assert(quo);
	assert(rem);
	if (denum->size == 0) {
		bigint_errno = BIGINT_ERR_DIV_BY_ZERO;
		ELOG_STR("DIVISION BY ZERO");
		return NULL;
	}
	size_t rem_size;
	size_t quo_size = div_recpr_cap(num->size, denum->size, recpr->size, precision, &rem_size);
	*quo = bigint_reserve(quo, quo_size, CLEAR_DATA);
	if (!*quo) return NULL;
	*rem = bigint_reserve(rem, rem_size, CLEAR_DATA);
	if (!*rem) return NULL;
	(*quo)->size = div_recpr(AS_SLICE(num), AS_SLICE(denum), AS_SLICE(recpr), precision, (*quo)->data, (*rem)->data, &rem_size);
	(*rem)->size = rem_size;
	return *quo;
}

BigInt bigint_div_recpr(ConstBigInt num, ConstBigInt denum, ConstBigInt recpr, size_t precision, BigInt* quo, BigInt* rem) {
	assert(num);
	assert(denum);
	assert(recpr);
	assert(quo);
	assert(rem);
	bool quo_sign = num->sign ^ denum->sign;
	bool rem_sign = num->sign;
	BigInt out = bigint_udiv_recpr(num, denum, recpr, precision, quo, rem);
	if (!out) return NULL;
	out->sign = quo_sign;
	(*rem)->sign = rem_sign;
	return out;
}

size_t bigint_decimal_width(ConstBigInt num) {
	return decimal_width(bigint_width(num));
}

