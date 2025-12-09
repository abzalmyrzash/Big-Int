#pragma once
#include "bigint.h"
#include "bigint_def.h"
#include "bigint_alias.h"
#include <assert.h>
#include <intrin.h>
#include <limits.h>
#include <stddef.h>
#include <string.h>
#include <inttypes.h>

typedef struct {
	void* buf;
	size_t bufsize;
} Context;

constexpr DataBlock BLOCK_MAX_VALUE = (DataBlock)~0;

constexpr uint8_t DATA_OFFSET = offsetof(struct bigint, data); // offset of data in bytes

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

static inline USmallInt ABS(SmallInt v) {
	if (v < 0) return ~(USmallInt)v + 1;
	return (USmallInt)v;
}

// returns the low part of a * b, high part in *high
static inline DataBlock block_mul(DataBlock a, DataBlock b, DataBlock* high) {
	BigMult mul = (BigMult)a * b;
	*high = mul >> BLOCK_WIDTH;
	return (DataBlock)mul;
}

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

static Slice split(Slice* slice, size_t size) {
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

static inline size_t set_usmall(DataBlock* data, USmallInt small) {
	if (small == 0) return 0;
	data[0] = small;
	return 1;
}

static inline size_t set_small(DataBlock* data, bool* sign, USmallInt small) {
	if (small == 0) return 0;
	*sign = small < 0;
	data[0] = ABS(small);
	return 1;
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
		if (!carry) {
			memmove(out + i, a.data + i, (a.size - i) * sizeof(DataBlock));
			return a.size;
		}
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

// out = a - b;
// a and b are unsigned but output can be negative;
static size_t sub(Slice a, Slice b, DataBlock* out, bool* out_sign) {
	int cmp_res = cmp(a, b);
	if (cmp_res == 0) {
		if (out_sign) *out_sign = 0;
		return 0;
	}
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
		size_t size = usub(a, b, out);
		assert(size == size_without_zeros(out, size));
		return size;
	} else {
		if (out_sign) *out_sign = 1;
		size_t size = usub(b, a, out);
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

static size_t sub_lshifted(Slice a, size_t lshift, Slice b, DataBlock* out, bool* out_sign) {
	int cmp_res = cmp_lshifted(a, lshift, b);
	if (cmp_res == 0) {
		if (out_sign) *out_sign = 0;
		return 0;
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
		size_t size = usub_lshifted(a, lshift, b, out);
		assert(size == size_without_zeros(out, size));
		return size;
	} else {
		if (out_sign) *out_sign = 1;
		size_t size = usub_b_lshifted(b, a, lshift, out);
		assert(size == size_without_zeros(out, size));
		return size;
	}
}

static size_t add_signed(Slice a, bool a_sign, Slice b, bool b_sign,
		DataBlock* out, bool* out_sign, bool assign_carry) {
	if (a_sign == b_sign) {
		*out_sign = a_sign;
		return add(a, b, out, assign_carry);
	}
	size_t size = sub(a, b, out, out_sign);
	*out_sign ^= a_sign;
	return size;
}

// (a << lshift) + b, signed
// NOTE: lshift by number of blocks not bits
static size_t add_lshifted_signed(Slice a, size_t lshift, bool a_sign,
		Slice b, bool b_sign, DataBlock* out, bool* out_sign, bool assign_carry) {
	if (a_sign == b_sign) {
		*out_sign = a_sign;
		return add_lshifted(a, lshift, b, out, assign_carry);
	}
	size_t size = sub_lshifted(a, lshift, b, out, out_sign);
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

// returns buffer size required for karatsuba in bytes
static inline size_t karatsuba_buffer_size(size_t n) {
	size_t sum = 0;
	static constexpr size_t threshold =
		BIGINT_KARATSUBA_THRESHOLD > 2 ? BIGINT_KARATSUBA_THRESHOLD : 2;
	while (n >= threshold) {
		if (n & 1) n += 1;
		sum += n;
		n >>= 1;
	}
	return sum * sizeof(DataBlock);
}

static int cnt_karatsuba = 0;
static int cnt_basicmul = 0;

static size_t mul(Slice a, Slice b, DataBlock* out, DataBlock* buffer) {
	const size_t n = a.size > b.size ? a.size : b.size;
	const size_t m = a.size < b.size ? a.size : b.size;
	size_t x;
	if (n >= BIGINT_KARATSUBA_THRESHOLD
			// standard algorithm is O(n * m)
			// karatsuba is O(n^1.57) or roughly O(n^1.5)
			// if m is greater than square root of n
			// then n * m > n * n^0.5 = n^1.5
			// so we will use m * m > n as another threshold for karatsuba
			// (also check for overflow)
			&& (__builtin_mul_overflow(m, m, &x) || x > n)) {
		cnt_karatsuba++;
		return karatsuba(a, b, n, out, buffer);
	} else {
		cnt_basicmul++;
		return mul_basic(a, b, out);
	}
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
// (a-b)(d-c)+ac+bd is used for the middle part
// because (a+b) and (c+d) can overflow
static size_t karatsuba(Slice b, Slice d, size_t n, DataBlock* out, DataBlock* buffer) {
	if(b.size == 0 || d.size == 0) return 0;
	assert(n == (b.size > d.size ? b.size : d.size));

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

	DataBlock* const abcd_data = out;
	DataBlock* const ab_data = buffer;
	DataBlock* const dc_data = buffer + m;
	DataBlock* const ac_data = buffer;
	DataBlock* const bd_data = buffer;
	DataBlock* const buf_data = buffer + 2 * m;

	bool ab_sign;
	bool dc_sign;
	bool abcd_sign;
	
	const Slice a_minus_b = {
		.data = ab_data,
		.size = sub(a, b, ab_data, &ab_sign)
	};
	assert(a_minus_b.size <= m);
	assert(a_minus_b.size == size_without_zeros(a_minus_b.data, a_minus_b.size));

	const Slice d_minus_c = {
		.data = dc_data,
		.size = sub(d, c, dc_data, &dc_sign)
	};
	assert(d_minus_c.size <= m);
	assert(d_minus_c.size == size_without_zeros(d_minus_c.data, d_minus_c.size));

	abcd_sign = ab_sign ^ dc_sign;

	// abcd = (a-b)(d-c)
	Slice abcd = {
		.data = abcd_data,
		.size = mul(a_minus_b, d_minus_c, abcd_data, buf_data)
	};
	assert(abcd.size <= 2 * m);
	assert(abcd.size == size_without_zeros(abcd.data, abcd.size));

	const Slice a_times_c = { 
		.data = ac_data,
		.size = mul(a, c, ac_data, buf_data)
	};
	assert(a_times_c.size <= 2 * (n - m));
	assert(a_times_c.size == size_without_zeros(a_times_c.data, a_times_c.size));

	// abcd = abcd + ac
	abcd.size = add_signed(abcd, abcd_sign, a_times_c, 0, abcd_data, &abcd_sign, true);
	assert(abcd.size == size_without_zeros(abcd.data, abcd.size));
	// abcd = (a_times_c << m) + abcd
	abcd.size = add_lshifted_signed(a_times_c, m, 0, abcd, abcd_sign, abcd_data, &abcd_sign, true);
	assert(abcd.size == size_without_zeros(abcd.data, abcd.size));

	const Slice b_times_d = {
		.data = bd_data,
		.size = mul(b, d, bd_data, buf_data)
	};
	assert(b_times_d.size <= 2 * m);
	assert(b_times_d.size == size_without_zeros(b_times_d.data, b_times_d.size));

	// abcd = abcd + bd
	abcd.size = add_signed(abcd, abcd_sign, b_times_d, 0, abcd_data, &abcd_sign, true);
	assert(abcd.size == size_without_zeros(abcd.data, abcd.size));
	// move abcd from out to out + m
	memmove(out + m, out, abcd.size * sizeof(DataBlock));
	abcd.data += m;
	// (abcd << m) + bd
	return add_lshifted_signed(abcd, m, abcd_sign, b_times_d, 0, out, &abcd_sign, true);
}

// out = a^pow % mod
static size_t pow_mod(Slice a, Slice pow, Slice mod, DataBlock* out, DataBlock* buffer) {
	const size_t out_cap = mod.size;
	Slice b = { .data = out, .size = 1 };
	out[0] = 1;
	const int bits = BLOCK_WIDTH - __builtin_clzg(pow.data[pow.size - 1]);
	for (int i = bits - 1; i >= 0; i--) {
		
	}
	for (size_t i = pow.size - 1; i > 0; i--) {
		DataBlock pow_block = pow.data[i - 1];
		for (int j = BLOCK_WIDTH - 1; j >= 0; j--) {
			// b = b * b
			b.size = mul(b, b, buffer, buffer + out_cap);
			memcpy(out, buffer, b.size * sizeof(DataBlock));
			// todo: take modulo
			if (pow_block & ((DataBlock)1 << j)) {
				// b = a * b
				b.size = mul(a, b, buffer, buffer + out_cap);
				memcpy(out, buffer, b.size * sizeof(DataBlock));
				// todo: take modulo
			}
		}
	}
	return b.size;
}

static size_t newton_reciprocal() {
	
}

static size_t copy(Slice z, DataBlock* out) {
	if (z.data != out) {
		memmove(out, z.data, z.size * sizeof(DataBlock));
	}
	return z.size;
}

static size_t lshift_blocks(Slice z, size_t shift, DataBlock* out) {
	if (z.size == 0) return 0;
	if (shift == 0) return copy(z, out);
	memset(out, 0, shift * sizeof(DataBlock));
	memmove(out + shift, z.data, z.size * sizeof(DataBlock));
	return z.size + shift;
}

static size_t rshift_blocks(Slice z, size_t shift, DataBlock* out) {
	if (z.size <= shift) return 0;
	memmove(out, z.data + shift, (z.size - shift) * sizeof(DataBlock));
	return z.size - shift;
}

static size_t lshift(Slice z, size_t shift, DataBlock* out) {
	if (z.size == 0) return 0;

	const size_t shift_blocks = shift / BLOCK_WIDTH;
	const size_t shift_bits = shift % BLOCK_WIDTH;
	if (shift_bits == 0) return lshift_blocks(z, shift_blocks, out);

	const DataBlock LOW_BITS = (1ULL << (BLOCK_WIDTH - shift_bits)) - 1;
	const DataBlock HIGH_BITS = ~0ULL - LOW_BITS;

	size_t out_size = z.size + shift_blocks;

	size_t i = z.size - 1;
	DataBlock hi, lo;

	hi = z.data[i] & HIGH_BITS;
	lo = z.data[i] & LOW_BITS;
	if (hi) {
		out_size++;
		out[i + shift_blocks + 1] = hi >> (BLOCK_WIDTH - shift_bits);
	}
	out[i + shift_blocks] = lo << shift_bits;

	if (i--) while (1) {
		hi = z.data[i] & HIGH_BITS;
		lo = z.data[i] & LOW_BITS;
		out[i + shift_blocks] = lo << shift_bits;
		out[i + shift_blocks + 1] |= hi >> (BLOCK_WIDTH - shift_bits);
		if (i == 0) break;
		i--;
	}

	memset(out, 0, shift_blocks * sizeof(DataBlock));
	return out_size;
}

static size_t rshift(Slice z, size_t shift, DataBlock* out) {
	size_t shift_blocks = shift / BLOCK_WIDTH;
	size_t shift_bits = shift % BLOCK_WIDTH;
	if (shift_bits == 0) return rshift_blocks(z, shift_blocks, out);

	if (shift_blocks >= z.size) return 0;

	size_t out_size = z.size - shift_blocks;
	if (shift_bits >= BLOCK_WIDTH - CLZ((DataBlock)z.data[z.size - 1])) {
		if (out_size <= 1) return 0;
		out_size--;
	}

	const DataBlock LOW_BITS = (1ULL << shift_bits) - 1;
	const DataBlock HIGH_BITS = ~0ULL - LOW_BITS;

	size_t i;
	DataBlock hi, lo;
	for (i = 0; i < out_size - 1; i++) {
		hi = z.data[i + shift_blocks] & HIGH_BITS;
		lo = z.data[i + shift_blocks + 1] & LOW_BITS;
		out[i] = (hi >> shift_bits) | lo << (BLOCK_WIDTH - shift_bits);
	}

	hi = z.data[i + shift_blocks] & HIGH_BITS;
	if (i + shift_blocks + 1 < z.size)
		lo = z.data[i + shift_blocks + 1] & LOW_BITS;
	else lo = 0;
	out[i] = (hi >> shift_bits) | lo << (BLOCK_WIDTH - shift_bits);

	return out_size;
}

static inline void bit_set(DataBlock* z, size_t bitpos) {
	z[bitpos / BLOCK_WIDTH] |= (1ULL << (bitpos % BLOCK_WIDTH));
}

static inline void bit_unset(DataBlock* z, size_t bitpos) {
	z[bitpos / BLOCK_WIDTH] &= ~(1ULL << (bitpos % BLOCK_WIDTH));
}

static inline void bit_toggle(DataBlock* z, size_t bitpos) {
	z[bitpos / BLOCK_WIDTH] ^= (1ULL << (bitpos % BLOCK_WIDTH));
}

static size_t div_basic(const Slice a, const Slice b, DataBlock* out, DataBlock* rem_data, CapField* rem_size, DataBlock* buffer) {
	const size_t size_diff = a.size - b.size;
	DataBlock* c_data = buffer;
	Slice rem = {
		.data = rem_data,
		.size = copy(a, rem_data)
	};
	Slice c = {
		.data = c_data,
		.size = lshift_blocks(b, size_diff, c_data)
	};

	int cmp_res = cmp(rem, c);
	if (cmp_res == 0) {
		*rem_size = 0;
		out[size_diff] = 1;
		memset(out, 0, size_diff * sizeof(DataBlock));
		return size_diff + 1;
	}
	if (cmp_res < 0 && size_diff == 0) {
		*rem_size = rem.size;
		return 0;
	}

	size_t bitpos;
	size_t out_size;
	const int c_clz = __builtin_clzg(c.data[c.size - 1]);
	const int rem_clz = __builtin_clzg(rem.data[rem.size - 1]);
	const int clz_diff = c_clz - rem_clz;

	if (cmp_res < 0) {
		assert(clz_diff <= 0);
		out_size = size_diff;
		memset(out, 0, out_size * sizeof(DataBlock));
		c.size = rshift(c, -clz_diff, c_data);
	} else {
		assert(clz_diff >= 0);
		out_size = size_diff + 1;
		memset(out, 0, out_size * sizeof(DataBlock));
		c.size = lshift(c, clz_diff, c_data);
	}
	bitpos = size_diff * BLOCK_WIDTH + clz_diff;

	cmp_res = cmp(rem, c);
	if (cmp_res < 0) {
		c.size = rshift(c, 1, c_data);
		bitpos--;
	}
	rem.size = usub(rem, c, rem_data);
	bit_set(out, bitpos);

	while(cmp(rem, b) >= 0) {
		while (cmp(rem, c) < 0) {
			c.size = rshift(c, 1, c_data);
			assert(bitpos > 0);
			bitpos--;
		}
		rem.size = usub(rem, c, rem_data);
		bit_set(out, bitpos);
	}
	*rem_size = rem.size;

	return out_size;
}
