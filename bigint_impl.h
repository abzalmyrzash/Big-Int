#pragma once
#include "bigint.h"
#include "bigint_def.h"
#include "bigint_alias.h"
#include "bigint_impl_basic.h"
#include <assert.h>
#include <immintrin.h>
#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <string.h>
#include <inttypes.h>

#define POSITIVE 0
#define NEGATIVE 1

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

#define MUT(x) _Generic((x), \
	const Block*: (Block*)(x),   \
	default: assert(0 && "MUT: invalid argument!\n"))

#define MUT_DATA(x) _Generic((x), \
	Slice: (Block*)(x).data, \
	default: assert(0 && "MUT_DATA: invalid argument!\n"))

typedef struct {
	Block* buf;
	size_t bufsize;
} Context;

static constexpr Block VALUE_TEN = 10;
static const Slice TEN = {
	.data = (Block*)&VALUE_TEN,
	.size = 1
};

constexpr Block BLOCK_MAX_VALUE = (Block)~0;

constexpr uint8_t DATA_OFFSET = offsetof(struct bigint, data); // offset of data in bytes

#define HEX_PREFIX "0x"

#define IS_SIGNED_TYPE(TYPE)    !IS_UNSIGNED_TYPE(TYPE)
#define IS_UNSIGNED_TYPE(TYPE)  ((TYPE)-1 > 0)
#define MAX_DEC_INT_DIGITS(TYPE)                \
  ((((sizeof(TYPE) * CHAR_BIT * 1233) >> 12) + 1) \
    + IS_SIGNED_TYPE(TYPE))

#define HEX_DIGITS_PER_BLOCK (sizeof(Block) * 2)
#define MAX_DEC_DIGITS_PER_BLOCK MAX_DEC_INT_DIGITS(Block)
#define DEC_DIGITS_PER_SPACE 3

constexpr int HEX_PREFIX_LEN = sizeof(HEX_PREFIX) - 1;

#define CLZ(x) __builtin_clzg(x)
#define CTZ(x) __builtin_ctzg(x)

#if BIGINT_BLOCK_WIDTH == 64
#ifdef _WIN32
	#define ADDCARRY _addcarry_u64
	#define SUBBORROW _subborrow_u64
#elif __linux__
	#define ADDCARRY(carry, a, b, out) _addcarry_u64(carry, a, b, (unsigned long long*)out)
	#define SUBBORROW(carry, a, b, out) _subborrow_u64(carry, a, b, (unsigned long long*)out)
#endif
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
static inline Block block_mul(Block a, Block b, Block* high) {
	BigMult mul = (BigMult)a * b;
	*high = mul >> BLOCK_WIDTH;
	return (Block)mul;
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

// compare (a <<< lshift) with b
// _sls stands for superleftshifted
static int cmp_sls(Slice a, size_t lshift, Slice b) {
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

static void hex_to_bin(Block x, char* str) {
	for (int i = 0; i < BLOCK_WIDTH; i++) {
		str[i] = (x & (1ULL << (BLOCK_WIDTH - i - 1))) ? '1' : '0';
	}
	str[BLOCK_WIDTH] = '\0';
}

static void print_slice_bin(Slice a) {
	if (a.size == 0) {
		printf("0");
		return;
	}
	static char bin[BLOCK_WIDTH + 1] = { 0 };
	hex_to_bin(a.data[a.size - 1], bin);
	printf("%s", bin);
	a.size--;
	while(a.size > 0) {
		hex_to_bin(a.data[a.size - 1], bin);
		printf("%s", bin);
		a.size--;
	}
}

static void insert_char(char* str, char ch, size_t pos) {
	int len = strlen(str);
	memmove(str + pos + 1, str + pos, len - pos);
	str[pos] = ch;
}

static void print_slice_bin_point(Slice a, size_t point) {
	if (a.size == 0) {
		printf("0");
		return;
	}
	size_t point_block = point / BLOCK_WIDTH;
	size_t point_pos = BLOCK_WIDTH - point % BLOCK_WIDTH;
	static char bin[BLOCK_WIDTH + 2] = { 0 };
	hex_to_bin(a.data[a.size - 1], bin);
	if (point_block == a.size - 1) {
		insert_char(bin, '.', point_pos);
	}
	printf("%s", bin);
	a.size--;
	while(a.size > 0) {
		hex_to_bin(a.data[a.size - 1], bin);
		if (point_block == a.size - 1) {
			insert_char(bin, '.', point_pos);
		}
		printf(" %s", bin);
		a.size--;
	}
}

// size in bits
static inline size_t width(Slice a) {
	if (a.size == 0) return 0;
	return a.size * BLOCK_WIDTH - CLZ(a.data[a.size - 1]);
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

static inline size_t set_usmall(Block* data, USmallInt small) {
	if (small == 0) return 0;
	data[0] = small;
	return 1;
}

static inline size_t set_small(Block* data, bool* sign, SmallInt small) {
	if (small == 0) return 0;
	*sign = small < 0;
	data[0] = ABS(small);
	return 1;
}

// a.size must be >= b.size
static size_t add_a(Slice a, Slice b, Block* out, bool assign_carry) {
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
			memmove(out + i, a.data + i, (a.size - i) * sizeof(Block));
			return a.size;
		}
		carry = ADDCARRY(carry, a.data[i], 0, &out[i]);
	}
	if (assign_carry && carry) out[a.size] = carry;
	return a.size + carry;
}

static size_t add(Slice a, Slice b, Block* out, bool assign_carry) {
	if (a.size >= b.size) {
		return add_a(a, b, out, assign_carry);
	}
	else {
		return add_a(b, a, out, assign_carry);
	}
}

// out = (a <<< lshift) + b;
static size_t add_sls(Slice a, size_t lshift, Slice b, Block* out, bool assign_carry) {
	if (a.size == 0) {
		memmove(out, b.data, b.size * sizeof(Block));
		return b.size;
	}
	if (b.size <= lshift) {
		memmove(out, b.data, b.size * sizeof(Block));
		memset(out + b.size, 0, (lshift - b.size) * sizeof(Block));
		memmove(out + lshift, a.data, a.size * sizeof(Block));
		return a.size + lshift;
	}
	memmove(out, b.data, lshift * sizeof(Block));
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
static size_t usub(Slice a, Slice b, Block* out) {
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
static size_t sub(Slice a, Slice b, Block* out, bool* out_sign) {
	int cmp_res = cmp(a, b);
	if (cmp_res == 0) {
		if (out_sign) *out_sign = 0;
		return 0;
	}
	if (a.size == 0) {
		*out_sign = 1;
		memmove(out, b.data, b.size * sizeof(Block));
		return b.size;
	}
	if (b.size == 0) {
		*out_sign = 0;
		memmove(out, a.data, a.size * sizeof(Block));
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

// out = (a <<< lshift) - b;
// Requirements: (a <<< lshift) >= b;
// _sls stands for superleftshifted
static size_t usub_sls(Slice a, size_t lshift, Slice b, Block* out) {
	assert(cmp_sls(a, lshift, b) >= 0);
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
// _b_sls stands for b is superleftshifted
static size_t usub_b_sls(Slice a, Slice b, size_t lshift, Block* out) {
	assert(cmp_sls(b, lshift, a) <= 0);
	memmove(out, a.data, lshift * sizeof(Block));
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

// (a <<< lshift) - b
// _sls stands for superleftshifted
static size_t sub_sls(Slice a, size_t lshift, Slice b, Block* out, bool* out_sign) {
	int cmp_res = cmp_sls(a, lshift, b);
	if (cmp_res == 0) {
		if (out_sign) *out_sign = 0;
		return 0;
	}
	if (a.size == 0) {
		*out_sign = 1;
		memmove(out, b.data, b.size * sizeof(Block));
		return b.size;
	}
	if (b.size == 0) {
		*out_sign = 0;
		memset(out, 0, lshift * sizeof(Block));
		memmove(out + lshift, a.data, a.size * sizeof(Block));
		return a.size + lshift;
	}
	if (cmp_res > 0) {
		if (out_sign) *out_sign = 0;
		size_t size = usub_sls(a, lshift, b, out);
		assert(size == size_without_zeros(out, size));
		return size;
	} else {
		if (out_sign) *out_sign = 1;
		size_t size = usub_b_sls(b, a, lshift, out);
		assert(size == size_without_zeros(out, size));
		return size;
	}
}

static size_t add_signed(Slice a, bool a_sign, Slice b, bool b_sign,
		Block* out, bool* out_sign, bool assign_carry) {
	if (a_sign == b_sign) {
		*out_sign = a_sign;
		return add(a, b, out, assign_carry);
	}
	size_t size = sub(a, b, out, out_sign);
	*out_sign ^= a_sign;
	return size;
}

// (a <<< lshift) + b, signed
// _sls stands for superleftshifted
static size_t add_sls_signed(Slice a, size_t lshift, bool a_sign,
		Slice b, bool b_sign, Block* out, bool* out_sign, bool assign_carry) {
	if (a_sign == b_sign) {
		*out_sign = a_sign;
		return add_sls(a, lshift, b, out, assign_carry);
	}
	size_t size = sub_sls(a, lshift, b, out, out_sign);
	*out_sign ^= a_sign;
	return size;
}

static void long_mul_add(Slice a, Slice b, Block* out) {
	if (a.size == 0 || b.size == 0) return;
	for (size_t i = 0; i < a.size; i++) {
		Block d1 = a.data[i];
		for (size_t j = 0; j < b.size; j++) {
			Block d2 = b.data[j];
			Block high;
			Block low = block_mul(d1, d2, &high);
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
}

static size_t long_mul(Slice a, Slice b, Block* out) {
	if (a.size == 0 || b.size == 0) return 0;
	const size_t max_size = a.size + b.size;
	memset(out, 0, max_size * sizeof(Block));
	for (size_t i = 0; i < a.size; i++) {
		Block d1 = a.data[i];
		for (size_t j = 0; j < b.size; j++) {
			Block d2 = b.data[j];
			Block high;
			Block low = block_mul(d1, d2, &high);
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

static size_t karatsuba(Slice b, Slice d, size_t n, Block* out, Block* buffer);

// returns buffer size required for multiplication
static inline size_t mul_buffer_size(size_t n) {
	size_t sum = 0;
	static constexpr size_t threshold =
		BIGINT_KARATSUBA_THRESHOLD > 2 ? BIGINT_KARATSUBA_THRESHOLD : 2;
	while (n >= threshold) {
		if (n & 1) n += 1;
		sum += n;
		n >>= 1;
	}
	return sum;
}

static int cnt_karatsuba = 0;
static int cnt_longmul = 0;

static size_t mul(Slice a, Slice b, Block* out, Block* buffer) {
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
		// cnt_karatsuba++;
		return karatsuba(a, b, n, out, buffer);
	} else {
		// cnt_longmul++;
		return long_mul(a, b, out);
	}
}

// Karatsuba multiplication
// @params:
// b - first number (or slice of a number);
// d - second number (or slice);
// n - max(b.size, d.size);
// out - filled with b * d, must be allocated outside
//       (at least (b.size + d.size) * sizeof(Block) bytes);
// buffer - memory for storing intermediate calculations
//          (karatsuba_buffer_size() * sizeof(Block) bytes)
// @return size of out
//
// Note: a slightly modified formula -
// (a-b)(d-c)+ac+bd is used for the middle part
// because (a+b) and (c+d) can overflow
static size_t karatsuba(Slice b, Slice d, size_t n, Block* out, Block* buffer) {

	if(b.size == 0 || d.size == 0) return 0;
	assert(n == (b.size > d.size ? b.size : d.size));

	if (n == 1) {
		Block high = 0;
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

	Block* const abcd_data = out;
	Block* const ab_data = buffer;
	Block* const dc_data = buffer + m;
	Block* const ac_data = buffer;
	Block* const bd_data = buffer;
	Block* const buf_data = buffer + 2 * m;

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
	abcd.size = add_sls_signed(a_times_c, m, 0, abcd, abcd_sign, abcd_data, &abcd_sign, true);
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
	memmove(out + m, out, abcd.size * sizeof(Block));
	abcd.data += m;
	// (abcd << m) + bd
	return add_sls_signed(abcd, m, abcd_sign, b_times_d, 0, out, &abcd_sign, true);
}

static size_t copy(Slice z, Block* out) {
	if (z.data != out) {
		memmove(out, z.data, z.size * sizeof(Block));
	}
	return z.size;
}

// super-left-shifts z, i.e. left-shifts by number of blocks rather than bits
// this operation is denoted as <<< (as opposed to << for left-shift by bits)
static size_t super_lshift(Slice z, size_t shift, Block* out) {
	if (z.size == 0) return 0;
	if (shift == 0) return copy(z, out);
	memmove(out + shift, z.data, z.size * sizeof(Block));
	memset(out, 0, shift * sizeof(Block));
	return z.size + shift;
}

// super-right-shifts z, i.e. right-shifts by number of blocks rather than bits
// this operation is denoted as >>> (as opposed to >> for right-shift by bits)
static size_t super_rshift(Slice z, size_t shift, Block* out) {
	if (z.size <= shift) return 0;
	memmove(out, z.data + shift, (z.size - shift) * sizeof(Block));
	return z.size - shift;
}

static size_t lshift(Slice z, size_t shift, Block* out) {
	if (z.size == 0) return 0;

	const size_t shift_blocks = shift / BLOCK_WIDTH;
	const int shift_bits = shift % BLOCK_WIDTH;
	if (shift_bits == 0) return super_lshift(z, shift_blocks, out);

	const Block LOW_BITS = (1ULL << (BLOCK_WIDTH - shift_bits)) - 1;
	const Block HIGH_BITS = ~0ULL - LOW_BITS;

	size_t out_size = z.size + shift_blocks;

	size_t i = z.size - 1;
	Block hi, lo;

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

	memset(out, 0, shift_blocks * sizeof(Block));
	return out_size;
}

static size_t rshift(Slice z, size_t shift, Block* out) {
	const size_t shift_blocks = shift / BLOCK_WIDTH;
	const int shift_bits = shift % BLOCK_WIDTH;
	if (shift_bits == 0) return super_rshift(z, shift_blocks, out);

	if (shift_blocks >= z.size) return 0;

	size_t out_size = z.size - shift_blocks;
	if (shift_bits >= BLOCK_WIDTH - CLZ((Block)z.data[z.size - 1])) {
		if (out_size <= 1) return 0;
		out_size--;
	}

	const Block LOW_BITS = (1ULL << shift_bits) - 1;
	const Block HIGH_BITS = ~0ULL - LOW_BITS;

	size_t i;
	Block hi, lo;
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

static inline void bit_set(Block* z, size_t bitpos) {
	z[bitpos / BLOCK_WIDTH] |= (1ULL << (bitpos % BLOCK_WIDTH));
}

static inline void bit_unset(Block* z, size_t bitpos) {
	z[bitpos / BLOCK_WIDTH] &= ~(1ULL << (bitpos % BLOCK_WIDTH));
}

static inline void bit_toggle(Block* z, size_t bitpos) {
	z[bitpos / BLOCK_WIDTH] ^= (1ULL << (bitpos % BLOCK_WIDTH));
}

// divides a / b, fills in the quotient in quo, remainder in rem_data and rem_size,
// uses buffer for temporary storage (up to (b.size + 1) * sizeof(Block) bytes)
// returns quotient's size
// uses basic long division algorithm
static size_t long_div(const Slice a, const Slice b, Block* quo,
		Block* rem_data, CapField* rem_size, Block* buffer) {
	assert(b.size > 0);

	if (a.size < b.size) {
		*rem_size = copy(a, rem_data);
		return 0;
	}

	const size_t size_diff = a.size - b.size;
	Block* c_data = buffer;

	// compare b <<< size_diff to a
	int cmp_res = cmp_sls(b, size_diff, a);

	// if b <<< size_diff == a, a / b = 1 <<< size_diff and a % b = 0
	if (cmp_res == 0) {
		*rem_size = 0;
		if (quo) {
			quo[size_diff] = 1;
			memset(quo, 0, size_diff * sizeof(Block));
		}
		return size_diff + 1;
	}

	// if (b <<< 0) = b > a, a / b = 0 and a % b = a
	if (cmp_res > 0 && size_diff == 0) {
		*rem_size = copy(a, rem_data);
		return 0;
	}

	// quotient size, we return this in the end
	size_t quo_size = 0;

	// this is a classic long division algorithm,
	// where we will repeatedly perform subtraction:
	// rem -= (b << bitpos)
	// and also set the bitpos-th bit in out
	// with the largest possible bitpos where rem >= (b << bitpos)
	size_t bitpos;
	// to be more efficient and save memory,
	// we store b leftshifted by up to BLOCK_WIDTH bits in c,
	// which leaves some integer number of whole blocks that weren't shifted
	// so we do subtraction by using usub_b_sls(rem, c, block_shift, rem_data)
	// which simulates subtracting c super-leftshifted by block_shift
	// without actually shifting c
	size_t block_shift = size_diff;
	// first we need to find the largest left shift we need to do
	// then subtract it from rem
	// and do the cycle of right shifting c and subtracting from rem
	// until rem becomes smaller than b

	// rem will be equal to a at first
	Slice rem = {
		.data = rem_data,
		.size = copy(a, rem_data)
	};

	// we will fill c_data by b shifted by some bits
	Slice c = {
		.data = c_data,
	};

	// number of leading zeros in a and b
	const int a_clz = CLZ(a.data[a.size - 1]);
	const int b_clz = CLZ(b.data[b.size - 1]);
	const int clz_diff = a_clz - b_clz;

	// first we shift b (into c) so that the clz of a and c match up

	// (b <<< size_diff) > a, rightshift by clz_diff
	if (cmp_res > 0) {
		assert(clz_diff >= 0 && clz_diff < BLOCK_WIDTH);
		if (quo) {
			quo_size = size_diff;
			memset(quo, 0, quo_size * sizeof(Block));
		}
		if (clz_diff > 0) {
			// rightshift b by clz_diff, however if we do it directly
			// we will lose some bits of data.
			// so instead we leftshift by BLOCK_WIDTH - clz_diff
			// and decrease block_shift
			// (note: clz_diff can only be >= 0 and
			//  if it's == 0 then we just copy b of course)
			c.size = lshift(b, BLOCK_WIDTH - clz_diff, c_data);
			block_shift--;
		}
		else {
			c.size = copy(b, c_data);
		}
	}

	// (b <<< size_diff) < a, leftshift by -clz_diff
	else {
		assert(clz_diff <= 0 && clz_diff > -(int)BLOCK_WIDTH);
		if (quo) {
			quo_size = size_diff + 1;
			memset(quo, 0, quo_size * sizeof(Block));
		}
		c.size = lshift(b, -clz_diff, c_data);
	}

	assert(CLZ(c.data[c.size - 1]) == a_clz);

	bitpos = size_diff * BLOCK_WIDTH - clz_diff;

	// (c <<< block_shift) > rem
	cmp_res = cmp_sls(c, block_shift, rem);
	if (cmp_res > 0) {
		// we can lose data if bitpos % BLOCK_WIDTH == 0
		if (bitpos % BLOCK_WIDTH == 0) {
			// rightshift by 1, but same trick as before
			c.size = lshift(c, BLOCK_WIDTH - 1, c_data);
			block_shift--;
		}
		else c.size = rshift(c, 1, c_data);
		bitpos--;
	}
	rem.size = usub_b_sls(rem, c, block_shift, rem_data);
	if (quo) bit_set(quo, bitpos);

	// subtract c from rem until rem becomes smaller than b
	while(cmp(rem, b) >= 0) {
		// right shift until c becomes smaller than rem
		do {
			// we can lose data if bitpos % BLOCK_WIDTH == 0
			if (bitpos % BLOCK_WIDTH == 0) {
				// rightshift by 1, but same trick as before
				c.size = lshift(c, BLOCK_WIDTH - 1, c_data);
				block_shift--;
			}
			else c.size = rshift(c, 1, c_data);
			assert(bitpos > 0);
			bitpos--;
		} while (cmp_sls(c, block_shift, rem) > 0);
		rem.size = usub_b_sls(rem, c, block_shift, rem_data);
		if (quo) bit_set(quo, bitpos);
	}
	*rem_size = rem.size;

	return quo_size;
}

static size_t count_trailing_zeros(Slice x) {
	for (size_t i = 0; i < x.size; i++) {
		if (x.data[i] != 0) {
			return i * BLOCK_WIDTH + CTZ(x.data[i]);
		}
	}
	return x.size * BLOCK_WIDTH;
}

static const Block VALUE_17_CONST = 17;
static const Slice _17_CONST = { .data = (Block*)&VALUE_17_CONST, .size = 1 };
static Slice  _48_OVER_17 = { NULL, 0 };
static size_t _48_OVER_17_CAP = 0;
static size_t _48_OVER_17_POINT = 0;
static size_t _48_OVER_17_PRECISION = 0;
static Slice  _48_OVER_17_REM = { NULL, 0 };
static size_t _48_OVER_17_REM_CAP = 0;
static Slice  _32_OVER_17 = { NULL, 0 };
static size_t _32_OVER_17_CAP = 0;
static size_t _32_OVER_17_POINT = 0;
static size_t _32_OVER_17_PRECISION = 0;
static Slice  _32_OVER_17_REM = { NULL, 0 };
static size_t _32_OVER_17_REM_CAP = 0;
static size_t NEWTON_PRECISION = 0;
static int NEWTON_STEPS = 0;

#define _17_WIDTH 5
#define _48_OVER_17_INT_WIDTH 2
#define _32_OVER_17_INT_WIDTH 1

static size_t newton_reciprocal_size(Slice d, size_t* buf_size) {
	const size_t d_point = width(d);
	const size_t x_point = d_point + _32_OVER_17_POINT;
	const size_t final_x_point = (x_point << NEWTON_STEPS) + (d_point << NEWTON_STEPS) - d_point;
	const size_t final_x_size = final_x_point / BLOCK_WIDTH + 1;

	const size_t prefinal_x_point = (x_point << (NEWTON_STEPS - 1)) + (d_point << (NEWTON_STEPS - 1)) - d_point;
	const size_t final_dx_size = (prefinal_x_point + d_point) / BLOCK_WIDTH + 1;

	const size_t new_x_cap = final_x_size + 1;
	const size_t dx_cap    = final_dx_size + 1;
	const size_t mul_cap   = mul_buffer_size(final_dx_size);

	*buf_size = new_x_cap + dx_cap + mul_cap;
	return final_x_size + 1;
}

// out = 1 / d
static size_t newton_reciprocal(Slice d, Block* out, Block* buffer) {
	assert(NEWTON_PRECISION > 0 && NEWTON_STEPS > 0);

	// we will "shift" d so it fits in range [0.5, 1)
	// by introducing a variable d_point representing the position of the radix point in d
	// it will be the width of d (size in bits)
	// e.g. 13 or binary 1101 has width 4 and will be "shifted" 4 bits,
	// becoming 13/(2^4)=0.8125 or binary 0.1101
	// similar variables are introduced for other numbers
	const size_t d_point = width(d);
	assert(d_point > 0);

	// calculate initial estimate of x = 48/17 - 32/17 * d
	Slice x = { .data = out };
	Block* tmp = buffer;
	// multiply 32/17 with d
	x.size = mul(_32_OVER_17, d, MUT_DATA(x), tmp);
	// x_point will be at d_point + 32/17_point when multiplying 32/17 * d
	size_t x_point = d_point + _32_OVER_17_POINT;
	// shift 48/17 so its point matches x_point
	Slice _4817_shf = { .data = buffer };
	_4817_shf.size = lshift(_48_OVER_17, x_point - _48_OVER_17_POINT, MUT_DATA(_4817_shf));
	print_slice_bin_point(_4817_shf, x_point);
	printf("\n");
	print_slice_bin_point(x, x_point);
	printf("\n");
	// x = 48/17 - 32/17 * d
	x.size = usub(_4817_shf, x, MUT_DATA(x));
	print_slice_bin_point(x, x_point);
	printf("\n");

	// x1 = 2 * x + d;
	// x2 = 2 * (2x + d) + d = 4x + 3d
	// x3 = 2 * (4x + d) + d = 8x + 7d
	// xn = 2^n * x + (2^n - 1) * d
	const size_t final_x_point = (x_point << NEWTON_STEPS) + (d_point << NEWTON_STEPS) - d_point;
	const size_t final_x_size = final_x_point / BLOCK_WIDTH + 1;
	printf("final_x_point = %zu\n", final_x_point);

	const size_t prefinal_x_point = (x_point << (NEWTON_STEPS - 1)) + (d_point << (NEWTON_STEPS - 1)) - d_point;
	const size_t last_mul_size = (prefinal_x_point + d_point) / BLOCK_WIDTH + 1;
	const size_t final_dx_size = (prefinal_x_point + d_point) / BLOCK_WIDTH + 1;

	const size_t new_x_cap = final_x_size + 1;
	const size_t dx_cap    = final_dx_size + 1;
	const size_t mul_cap   = mul_buffer_size(final_dx_size);
	printf("new_x_cap = %zu\n", new_x_cap);
	printf("dx_cap = %zu\n", dx_cap);
	printf("mul_cap = %zu\n", mul_cap);

	Slice new_x  = { .data = buffer };
	Slice dx     = { .data = new_x.data + new_x_cap };
	Block* mul_buf = MUT(dx.data + dx_cap);
	Block _1_data = 1;
	Slice _1 = { .data = &_1_data, .size = 1 };
	bool dx_sign;
	bool out_sign;
	printf("NEWTON_STEPS = %d\n", NEWTON_STEPS);

	for (int i = 0; i < NEWTON_STEPS; i++) {
		// x + x * (1 - d * x)
		// d * x
		dx.size = mul(d, x, MUT_DATA(dx), mul_buf);
		printf("dx =           ");
		print_slice_bin_point(dx, d_point + x_point);
		printf("\n");

		// shift 1 by d_point + x_point
		size_t shf = d_point + x_point;
		size_t shf_blocks = shf / BLOCK_WIDTH;
		size_t shf_bits = shf % BLOCK_WIDTH;
		_1_data = 1ULL << shf_bits;

		printf("%zu\n", shf);
		printf("_1 =           ");
		print_slice_bin_point(_1, d_point + x_point);
		printf("\n");

		// 1 - dx
		dx.size = sub_sls(_1, shf_blocks, dx, MUT_DATA(dx), &dx_sign);

		printf("1 - dx =       ");
		print_slice_bin_point(dx, d_point + x_point);
		printf("(%u)\n", dx.size);
		print_slice_bin_point(x, x_point);
		printf("(%u)\n", x.size);

		// x * (1 - dx)
		new_x.size = mul(x, dx, MUT_DATA(new_x), mul_buf);

		printf("x * (1 - dx) = ");
		print_slice_bin_point(new_x, d_point + 2 * x_point);
		printf("(%u)\n", new_x.size);

		printf("x = ");
		print_slice_bin_point(x, x_point);
		printf("(%u)\n", x.size);
		printf("d_point + x_point = %zu\n", d_point + x_point);
		x.size = lshift(x, d_point + x_point, MUT_DATA(x));
		printf("x = ");
		print_slice_bin_point(x, d_point + 2 * x_point);
		printf("(%u)\n", x.size);
		x.size = add_signed(x, POSITIVE, new_x, dx_sign, MUT_DATA(x), &out_sign, true);
		x_point += d_point + x_point;

		printf("x = ");
		print_slice_bin_point(x, x_point);
		printf("(%u)\n", x.size);
	}

	printf("\n");
	printf("\n");
	print_slice(x);
	return x.size;
}

static inline size_t powmod_buffer_size(size_t mod_size) {
	const size_t mul_cap = mod_size * 2 + mul_buffer_size(mod_size);
	const size_t div_cap = mod_size * 3 + 1;
	return mul_cap > div_cap ? mul_cap : div_cap;
}

// out = a^pow % mod
static size_t powmod(Slice a, Slice pow, Slice mod, Block* out_data, Block* buffer) {
	assert(cmp(a, mod) < 0);
	const size_t mul_cap = mod.size * 2;
	Block* const mul_buf = buffer + mul_cap;
	Block* const div_buf = buffer + mul_cap;
	Slice out = { .data = out_data, .size = 1 };
	Slice buf = { .data = buffer };
	out_data[0] = 1;

	for (size_t i = pow.size; i > 0; i--) {
		Block pow_block = pow.data[i - 1];
		const int bits = (i == pow.size) ?
			BLOCK_WIDTH - CLZ(pow_block) : BLOCK_WIDTH;
		for (int j = bits - 1; j >= 0; j--) {
			// buffer = out * out
			buf.size = mul(out, out, buffer, mul_buf);
			// out = buffer % mod
			long_div(buf, mod, NULL, buffer, &out.size, div_buf);
			memcpy(out_data, buffer, sizeof(Block) * out.size);
			if (pow_block & ((Block)1 << j)) {
				// buffer = a * out
				buf.size = mul(a, out, buffer, buffer + mul_cap);
				// out = buffer % mod
				long_div(buf, mod, NULL, buffer, &out.size, div_buf);
				memcpy(out_data, buffer, sizeof(Block) * out.size);
			}
		}
	}

	return out.size;
}

