#include "bigint_impl.h"
#include <assert.h>
#include <stdlib.h>
#include <stdarg.h>
#include "elog.h"
#include <errno.h>
#include "bigint_rand.h"

int bigint_errno = 0;

static thread_local Context ctx;

static unsigned long long cnt_alloc = 0;
static unsigned long long cnt_realloc = 0;
static unsigned long long cnt_free = 0;

static Block* ctx_buf_reserve(size_t cap);

#define AS_SLICE(x) (Slice) _Generic((x), \
	BigInt: (Slice){ .data = (x)->data, .size = (x)->size })

void bigint_init(void) {
	ctx.buf = NULL;
	ctx.bufsize = 0;

	INITIAL_NEWTON_PRECISION = log2(17);
	NEWTON_PRECISION = 4;
	NEWTON_STEPS = 0;

	_48_OVER_17.data = malloc(sizeof(Block));
	_32_OVER_17.data = malloc(sizeof(Block));
	_48_OVER_17_REM.data = malloc(sizeof(Block));
	_32_OVER_17_REM.data = malloc(sizeof(Block));

	_48_OVER_17_POINT = NEWTON_PRECISION - _48_OVER_17_INT_WIDTH;
	MUT_DATA(_48_OVER_17)[0] = (48ULL << _48_OVER_17_POINT) / 17;
	_48_OVER_17.size = 1;
	MUT_DATA(_48_OVER_17_REM)[0] = (48ULL << _48_OVER_17_POINT) % 17;
	_48_OVER_17_REM.size = 1;

	_32_OVER_17_POINT = NEWTON_PRECISION - _32_OVER_17_INT_WIDTH;
	MUT_DATA(_32_OVER_17)[0] = (32ULL << _32_OVER_17_POINT) / 17;
	_32_OVER_17.size = 1;
	MUT_DATA(_32_OVER_17_REM)[0] = (32ULL << _32_OVER_17_POINT) % 17;
	_32_OVER_17_REM.size = 1;

	size_t precision = 1000;
	BigInt rand_num = NULL;
	bigint_rand(&rand_num, precision);
	bigint_printf("%d\n", rand_num);

	set_newton_precision(precision);

	size_t bufsize;
	size_t size = newton_reciprocal_size(TEN, &bufsize);
	printf("size: %zu\n", size);
	printf("bufsize: %zu\n", bufsize);
	Block* buf = ctx_buf_reserve(bufsize);
	Slice recipr = { .data = malloc(sizeof(Block) * size) };
	recipr.size = newton_reciprocal(TEN, MUT_DATA(recipr), buf);
	print_slice(recipr);
	printf("(%u)\n", recipr.size);
	size_t mul_bufsize = mul_buffer_size(MAX(recipr.size, rand_num->size));
	buf = ctx_buf_reserve(mul_bufsize);
	BigInt mul_res = bigint_alloc(rand_num->size + recipr.size);
	mul_res->size = mul(AS_SLICE(rand_num), recipr, mul_res->data, buf);
	mul_res->size = rshift(AS_SLICE(mul_res), precision + width(TEN), mul_res->data);
	bigint_printf("%d\n", mul_res);

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
	free(ctx.buf);
	ctx.buf = NULL;
	ctx.bufsize = 0;
}

static size_t cnt_ctx_buf_realloc = 0;

static Block* ctx_buf_reserve(size_t cap) {
	if(ctx.bufsize < cap) {
		free(ctx.buf);
		ctx.buf = malloc(cap * sizeof(Block));
		if (!ctx.buf) {
			ELOG_CODE(errno);
			ELOG("Parameters: cap = %zu\n", cap);
			return NULL;
		}
		cnt_ctx_buf_realloc++;
		ctx.bufsize = cap;
	}
	return ctx.buf;
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
	printf("Buffer reallocations: %llu\n", cnt_ctx_buf_realloc);
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
	BigInt out = bigint_reserve(out_ptr, a->size, CLEAR_DATA);
	if (!out) return NULL;
	bool sign;
	out->size = sub(A, B, out->data, &sign);
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
	const size_t mul_buf_size = mul_buffer_size(n);

	Slice A = { .data = a->data, .size = a->size };
	Slice B = { .data = b->data, .size = b->size };

	Block* buf = ctx_buf_reserve(
			mul_buf_size
			+ (out == a) * A.size
			+ (out == b) * B.size);
	if (!buf) goto error;

	if (out == a) {
		A.data = memcpy(ctx.buf, A.data, A.size * sizeof(Block));
		buf += A.size;
	}
	else if (out == b) {
		B.data = memcpy(ctx.buf, B.data, B.size * sizeof(Block));
		buf += B.size;
	}

	out = bigint_reserve(&out, cap, CLEAR_DATA);
	if (!out) goto error;

	out->size = mul(A, B, out->data, buf);
	assert(out->size <= cap);
	assert(out->size == size_without_zeros(out->data, out->size));
	assert(({
		Block* check = malloc(cap * sizeof(Block));
		size_t check_size = long_mul(A, B, check);
		int cmp_res = cmp((Slice){out->data, out->size}, (Slice){check, check_size});
		free(check);
		cmp_res == 0;
	}));

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

	BigInt out = *out_ptr;
	BigInt rem = rem_ptr ? *rem_ptr : NULL;

	Slice A = { .data = a->data, .size = a->size };
	Slice B = { .data = b->data, .size = b->size };

	if (bigint_ucmp_small(b, 1) == 0) {
		out = bigint_copy(&out, a);
		if (!out) goto error;
		if (rem_ptr) {
			rem = bigint_set_zero(&rem);
			if (!rem) goto error;
		}
		goto ret;
	}

	int cmp_res = cmp(A, B);

	if (cmp_res < 0) {
		if (rem_ptr) {
			rem = bigint_copy(&rem, a);
			if (!rem) goto error;
		}
		out = bigint_set_zero(&out);
		if (!out) goto error;
		goto ret;
	}

	else if (cmp_res == 0) {
		if (rem_ptr) {
			rem = bigint_set_zero(&rem);
			if (!rem) goto error;
		}
		out = bigint_set_usmall(&out, 1);
		if (!out) goto error;
		goto ret;
	}
	
	size_t bufsize_blocks = b->size + 1;
	if (out == a) {
		bufsize_blocks += a->size;
	} else if (out == b) {
		bufsize_blocks += b->size;
	}
	if (rem == b) {
		bufsize_blocks += b->size;
	}
	else if (!rem_ptr) {
		bufsize_blocks += a->size;
	}

	Block* buf = ctx_buf_reserve(bufsize_blocks);
	if (!buf) {
		goto error;
	}

	Block* rem_data;
	CapField rem_size;

	if (out == a) {
		A.data = memcpy(buf, A.data, A.size * sizeof(Block));
		buf += a->size;
	} else if (out == b) {
		B.data = memcpy(buf, B.data, B.size * sizeof(Block));
		buf += b->size;
	}

	out = bigint_reserve(&out, A.size - B.size + 1, CLEAR_DATA);
	if (!out) goto error;

	if (rem_ptr) {
		rem = bigint_reserve(&rem, A.size, CLEAR_DATA);
		if (!rem) goto error;
		rem_data = rem->data;
		if (rem == b) {
			B.data = memcpy(buf, B.data, B.size * sizeof(Block));
			buf += b->size;
		}
	} else {
		rem_data = buf;
		buf += a->size;
	}

	out->size = long_div(A, B, out->data, rem_data, &rem_size, buf);
	if (rem) rem->size = rem_size;

ret:
	if (rem_ptr) *rem_ptr = rem;
	return *out_ptr = out;

error:
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

BigInt bigint_powmod(ConstBigInt a, ConstBigInt pow, ConstBigInt mod, BigInt* out_ptr) {
	assert(a);
	assert(pow);
	assert(mod);
	assert(out_ptr);
	
	BigInt out = bigint_reserve(out_ptr, mod->size, CLEAR_DATA);
	if (!out) return NULL;

	const size_t buf_size = powmod_buffer_size(mod->size);
	Block* buf = ctx_buf_reserve(buf_size);
	if (!buf) return NULL;

	Slice A = { a->data, a->size };
	Slice Pow = { pow->data, pow->size };
	Slice Mod = { mod->data, mod->size };
	out->size = powmod(A, Pow, Mod, out->data, buf);
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

	Slice Z = { z->data, z->size };

	const size_t q_cap = z->size;
	const size_t tmp_cap = z->size;
	static const size_t div_buf_cap = 2;

	const size_t bufsize = q_cap + tmp_cap + div_buf_cap;
	Block* const buf = ctx_buf_reserve(bufsize);
	if (!buf) return NULL;

	Block* const q_data = buf;
	Block* const tmp_data = q_data + q_cap;
	Block* const div_buf = tmp_data + tmp_cap;
	Slice tmp = {
		.data = tmp_data,
		.size = copy(Z, tmp_data)
	};

	USmallInt rem;
	CapField rem_size;

	while (tmp.size > 0) {
		tmp.size = long_div(tmp, TEN, q_data, tmp_data, &rem_size, div_buf);
		assert((rem_size == 1 && tmp_data[0] < 10) || rem_size == 0);
		rem = rem_size ? tmp_data[0] : 0;
		assert(cnt < max_digits + max_spaces || ({
					printf("%*s\n", (int)cnt, digits);
					0;
				}));
		digits[cnt++] = rem + '0';
		if (bifs.add_spaces &&
				(cnt + 1) % (DEC_DIGITS_PER_SPACE + 1) == 0) {
			assert(cnt < max_digits + max_spaces || ({
						printf("%*s\n", (int)cnt, digits);
						0;
					}));
			digits[cnt++] = ' ';
		}
		memcpy(tmp_data, q_data, tmp.size * sizeof(Block));
	}
	if (digits[cnt - 1] == ' ') cnt--;

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
	const size_t buf_size = pow10_cap + digit_cap + tmp_cap;
	Block* const buf = ctx_buf_reserve(buf_size);
	if (!buf) return NULL;

	Block* const pow10_data = buf;
	Block* const digit_data = pow10_data + pow10_cap;
	Block* const tmp_data = digit_data + digit_cap;
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

		digit.size = long_mul(tmp, pow10, digit_data);
		out->size = add((Slice){out->data, out->size}, digit, out->data, true);
		if (i == 0) break;
		i--;
		tmp.size = copy(pow10, tmp_data);
		pow10.size = long_mul(tmp, TEN, pow10_data);
	}

	return out;
}

bool set_newton_precision(size_t p) {
	NEWTON_STEPS = ceil( log2( (p + 1) / log2(17) ) );

	if (NEWTON_PRECISION >= p) return true;

	const size_t blocks = (p - 1) / BLOCK_WIDTH + 1;

	size_t dp = p - NEWTON_PRECISION;
	size_t dp_blocks = dp / BLOCK_WIDTH;
	size_t dp_bits = dp % BLOCK_WIDTH;
	size_t rem_cap = (dp + _17_WIDTH - 1) / BLOCK_WIDTH + 1;

	Block* tmp;
	if (_48_OVER_17_CAP < blocks) {
		tmp = realloc(MUT_DATA(_48_OVER_17), sizeof(Block) * blocks);
		if (!tmp) return false;
		_48_OVER_17.data = tmp;
		_48_OVER_17_CAP = blocks;
	}
	if (_32_OVER_17_CAP < blocks) {
		tmp = realloc(MUT_DATA(_32_OVER_17), sizeof(Block) * blocks);
		if (!tmp) return false;
		_32_OVER_17.data = tmp;
		_32_OVER_17_CAP = blocks;
	}

	if (_48_OVER_17_REM_CAP < rem_cap) {
		tmp = realloc(MUT_DATA(_48_OVER_17_REM), sizeof(Block) * rem_cap);
		if (!tmp) return false;
		_48_OVER_17_REM.data = tmp;
		_48_OVER_17_REM_CAP = rem_cap;
	}
	if (_32_OVER_17_REM_CAP < rem_cap) {
		tmp = realloc(MUT_DATA(_32_OVER_17_REM), sizeof(Block) * rem_cap);
		if (!tmp) return false;
		_32_OVER_17_REM.data = tmp;
		_32_OVER_17_REM_CAP = rem_cap;
	}

	_48_OVER_17.size = lshift(_48_OVER_17, dp, MUT_DATA(_48_OVER_17));
	Block _4817_mid_block = _48_OVER_17.data[dp_blocks];
	_32_OVER_17.size = lshift(_32_OVER_17, dp, MUT_DATA(_32_OVER_17));
	Block _3217_mid_block = _32_OVER_17.data[dp_blocks];

	Block* div_buf = ctx_buf_reserve(2);

	_48_OVER_17_REM.size = lshift(_48_OVER_17_REM, dp, MUT_DATA(_48_OVER_17_REM));
	long_div(_48_OVER_17_REM, _17_CONST,
			(Block*)_48_OVER_17.data, (Block*)_48_OVER_17_REM.data, &_48_OVER_17_REM.size, div_buf);
	MUT_DATA(_48_OVER_17)[dp_blocks] |= _4817_mid_block;
	_48_OVER_17_POINT += dp;

	_32_OVER_17_REM.size = lshift(_32_OVER_17_REM, dp, MUT_DATA(_32_OVER_17_REM));
	long_div(_32_OVER_17_REM, _17_CONST,
			(Block*)_32_OVER_17.data, (Block*)_32_OVER_17_REM.data, &_32_OVER_17_REM.size, div_buf);
	MUT_DATA(_32_OVER_17)[dp_blocks] |= _3217_mid_block;
	_32_OVER_17_POINT += dp;

	/*
	print_slice_bin_point(_48_OVER_17, _48_OVER_17_POINT);
	printf(":");
	print_slice(_48_OVER_17_REM);
	printf("\n");
	print_slice_bin_point(_32_OVER_17, _32_OVER_17_POINT);
	printf(":");
	print_slice(_32_OVER_17_REM);
	printf("\n");
	*/

	NEWTON_PRECISION = p;
	return true;
}

