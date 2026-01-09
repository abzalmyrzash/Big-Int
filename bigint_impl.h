#pragma once
#include "arena.h"
#include "bigint.h"
#include "bigint_def.h"
#include "bigint_alias.h"
#include "bigint_impl_basic.h"
#include "elog.h"
#include <assert.h>
#include <immintrin.h>
#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <string.h>
#include <inttypes.h>

#define POSITIVE 0
#define NEGATIVE 1

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

thread_local mem_arena* arena;

static constexpr Block VALUE_TEN = 10;
static const Slice TEN = {
	.data = (Block*)&VALUE_TEN,
	.size = 1
};

static constexpr Block VALUE_ONE = 1;
static const Slice ONE = {
	.data = (Block*)&VALUE_ONE,
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
	for (size_t i = 0; i < BLOCK_WIDTH; i++) {
		str[i] = (x & (1ULL << (BLOCK_WIDTH - i - 1))) ? '1' : '0';
	}
	str[BLOCK_WIDTH] = '\0';
}

static void print_slice_bin(Slice a) {
	if (a.size == 0) {
		printf("0");
		return;
	}
	char bin[BLOCK_WIDTH + 1] = { 0 };
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
	char bin[BLOCK_WIDTH + 3] = { 0 };
	Block first_block = a.data[a.size - 1];
	hex_to_bin(first_block, bin);
	if (point_block == a.size - 1) {
		insert_char(bin, '.', point_pos);
		if (a.size == 1) {
			char* str = bin + CLZ(first_block);
			if (str[0] == '.') insert_char(str, '0', 0);
			printf("%s", str);
			return;
		}
	}
	printf("%s", bin + CLZ(first_block));
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
		if (out_sign) *out_sign = POSITIVE;
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
		if (out_sign) *out_sign = POSITIVE;
		size_t size = usub(a, b, out);
		assert(size == size_without_zeros(out, size));
		return size;
	} else {
		if (out_sign) *out_sign = NEGATIVE;
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
		if (out_sign) *out_sign = POSITIVE;
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
		if (out_sign) *out_sign = POSITIVE;
		size_t size = usub_sls(a, lshift, b, out);
		assert(size == size_without_zeros(out, size));
		return size;
	} else {
		if (out_sign) *out_sign = NEGATIVE;
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

static size_t karatsuba(Slice b, Slice d, size_t n, Block* out);

static int cnt_karatsuba = 0;
static int cnt_longmul = 0;

static size_t mul(Slice a, Slice b, Block* out) {
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
		return karatsuba(a, b, n, out);
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
static size_t karatsuba(Slice b, Slice d, size_t n, Block* out) {

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

	bool ab_sign;
	bool dc_sign;
	bool abcd_sign;
	
	u64 restore_pos = arena->pos;
	Block* const ab_data = PUSH_ARRAY(arena, Block, m);
	const Slice a_minus_b = {
		.data = ab_data,
		.size = sub(a, b, ab_data, &ab_sign)
	};
	assert(a_minus_b.size <= m);
	assert(a_minus_b.size == size_without_zeros(a_minus_b.data, a_minus_b.size));

	Block* const dc_data = PUSH_ARRAY(arena, Block, m);
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
		.size = mul(a_minus_b, d_minus_c, abcd_data)
	};
	assert(abcd.size <= 2 * m);
	assert(abcd.size == size_without_zeros(abcd.data, abcd.size));
	arena_pop_to(arena, restore_pos);

	Block* const ac_data = PUSH_ARRAY(arena, Block, 2 * m);
	const Slice a_times_c = { 
		.data = ac_data,
		.size = mul(a, c, ac_data)
	};
	assert(a_times_c.size <= 2 * (n - m));
	assert(a_times_c.size == size_without_zeros(a_times_c.data, a_times_c.size));

	// abcd = abcd + ac
	abcd.size = add_signed(abcd, abcd_sign, a_times_c, 0, abcd_data, &abcd_sign, true);
	assert(abcd.size == size_without_zeros(abcd.data, abcd.size));
	// abcd = (a_times_c << m) + abcd
	abcd.size = add_sls_signed(a_times_c, m, 0, abcd, abcd_sign, abcd_data, &abcd_sign, true);
	assert(abcd.size == size_without_zeros(abcd.data, abcd.size));

	arena_pop_to(arena, restore_pos);

	Block* const bd_data = PUSH_ARRAY(arena, Block, 2 * m);
	const Slice b_times_d = {
		.data = bd_data,
		.size = mul(b, d, bd_data)
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
	size_t size = add_sls_signed(abcd, m, abcd_sign, b_times_d, 0, out, &abcd_sign, true);
	arena_pop_to(arena, restore_pos);
	return size;
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

static size_t round_rshift(Slice x, size_t shift, Block* out_data) {
	if (shift == 0) {
		memmove(out_data, x.data, sizeof(Block) * x.size);
		return x.size;
	}
	if (shift > width(x)) {
		return 0;
	}
	const size_t block = (shift - 1) / BLOCK_WIDTH;
	const size_t bit = (shift - 1) % BLOCK_WIDTH;
	const Block next_bit = out_data[block] & (1ULL << bit);
	const size_t size = rshift(x, shift, out_data);
	if (next_bit) {
		Slice out = { out_data, size };
		return add(out, ONE, out_data, true);
	}
	return size;
}

static size_t ceil_rshift(Slice x, size_t shift, Block* out_data) {
	if (shift == 0) {
		memmove(out_data, x.data, sizeof(Block) * x.size);
		return x.size;
	}
	if (shift > width(x)) {
		return 0;
	}
	const size_t block = shift / BLOCK_WIDTH;
	const size_t bit = shift % BLOCK_WIDTH;
	bool is_zero = bit ? x.data[block] & ((1ULL << bit) - 1) : true;
	if (is_zero) {
		for (size_t i = 0; i < block; i++) {
			if (x.data[i]) {
				is_zero = false;
				break;
			}
		}
	}
	const size_t size = rshift(x, shift, out_data);
	if (!is_zero) {
		Slice out = { out_data, size };
		return add(out, ONE, out_data, true);
	}
	return size;
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
// returns quotient's size
// uses basic long division algorithm
static size_t long_div(const Slice a, const Slice b, Block* quo,
		Block* rem_data, size_t* rem_size) {
	assert(b.size > 0);

	if (a.size < b.size) {
		*rem_size = copy(a, rem_data);
		return 0;
	}

	const size_t size_diff = a.size - b.size;

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

	size_t restore_pos = arena->pos;

	// we will fill c_data by b shifted by some bits
	Slice c = {
		.data = PUSH_ARRAY(arena, Block, b.size + 1)
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
			c.size = lshift(b, BLOCK_WIDTH - clz_diff, MUT_DATA(c));
			block_shift--;
		}
		else {
			c.size = copy(b, MUT_DATA(c));
		}
	}

	// (b <<< size_diff) < a, leftshift by -clz_diff
	else {
		assert(clz_diff <= 0 && clz_diff > -(int)BLOCK_WIDTH);
		if (quo) {
			quo_size = size_diff + 1;
			memset(quo, 0, quo_size * sizeof(Block));
		}
		c.size = lshift(b, -clz_diff, MUT_DATA(c));
	}

	assert(CLZ(c.data[c.size - 1]) == a_clz);

	bitpos = size_diff * BLOCK_WIDTH - clz_diff;

	// (c <<< block_shift) > rem
	cmp_res = cmp_sls(c, block_shift, rem);
	if (cmp_res > 0) {
		// we can lose data if bitpos % BLOCK_WIDTH == 0
		if (bitpos % BLOCK_WIDTH == 0) {
			// rightshift by 1, but same trick as before
			c.size = lshift(c, BLOCK_WIDTH - 1, MUT_DATA(c));
			block_shift--;
		}
		else c.size = rshift(c, 1, MUT_DATA(c));
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
				c.size = lshift(c, BLOCK_WIDTH - 1, MUT_DATA(c));
				block_shift--;
			}
			else c.size = rshift(c, 1, MUT_DATA(c));
			assert(bitpos > 0);
			bitpos--;
		} while (cmp_sls(c, block_shift, rem) > 0);
		rem.size = usub_b_sls(rem, c, block_shift, rem_data);
		if (quo) bit_set(quo, bitpos);
	}
	*rem_size = rem.size;

	arena_pop_to(arena, restore_pos);
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
static size_t _48_OVER_17_WIDTH = 0;
static Slice  _48_OVER_17_REM = { NULL, 0 };
static size_t _48_OVER_17_REM_CAP = 0;
static Slice  _32_OVER_17 = { NULL, 0 };
static size_t _32_OVER_17_CAP = 0;
static size_t _32_OVER_17_POINT = 0;
static size_t _32_OVER_17_WIDTH = 0;
static Slice  _32_OVER_17_REM = { NULL, 0 };
static size_t _32_OVER_17_REM_CAP = 0;
static double INITIAL_NEWTON_PRECISION = 0;
static int NEWTON_STEPS = 0;

static constexpr size_t _17_WIDTH = 5;
static constexpr size_t _48_OVER_17_INT_WIDTH = 2;
static constexpr size_t _32_OVER_17_INT_WIDTH = 1;

static size_t newton_reciprocal_cap(size_t precision) {
	return precision / BLOCK_WIDTH + 2;
}

// out = 1 / d
static size_t newton_reciprocal(Slice d, size_t precision, Block* out) {
	assert(precision > 0);

	const size_t x_cap = newton_reciprocal_cap(precision);

	// we will "shift" d so it fits in range [0.5, 1)
	// by introducing a variable d_point representing the position of the radix point in d
	// it will be the width of d (size in bits)
	// e.g. 13 or binary 1101 has width 4 and will be "shifted" 4 bits,
	// becoming 13/(2^4)=0.8125 or binary 0.1101
	// similar variables are introduced for other numbers
	const size_t d_point = width(d);
	assert(d_point > 0);
	assert(precision >= d_point);

	double p = INITIAL_NEWTON_PRECISION;
	size_t p_ceil = MIN(ceil(p), precision);

	u64 restore_pos = arena->pos;

	// calculate initial estimate of x = 48/17 - 32/17 * d
	Slice x = { .data = out };
	// multiply 32/17 with d
	x.size = mul(_32_OVER_17, d, MUT_DATA(x));
	// x_point will be at d_point + 32/17_point when multiplying 32/17 * d
	size_t x_point = d_point + _32_OVER_17_POINT;
	// shift 48/17 so its point matches x_point
	Slice _4817_shf = { .data = PUSH_ARRAY(arena, Block, x_point / BLOCK_WIDTH + 2) };
	_4817_shf.size = lshift(_48_OVER_17, x_point - _48_OVER_17_POINT, MUT_DATA(_4817_shf));
	// x = 48/17 - 32/17 * d
	
	x.size = usub(_4817_shf, x, MUT_DATA(x));

	arena_pop_to(arena, restore_pos);

	/*
	print_slice_bin_point(x, x_point);
	printf("\n");
	*/
	if (x_point > p_ceil) {
		x.size = rshift(x, x_point - p_ceil, MUT_DATA(x));
		x_point = p_ceil;
	}

	/*
	print_slice_bin_point(x, x_point);
	printf("\n");
	*/

	const size_t dx_cap  = x_cap + d.size;
	const size_t xdx_cap = x_cap + dx_cap;

	bool dx_sign;
	Slice dx       = { .data = PUSH_ARRAY(arena, Block, dx_cap) };
	Slice xdx      = { .data = PUSH_ARRAY(arena, Block, xdx_cap) };

	Block _1_data = 1;
	Slice _1 = { .data = &_1_data, .size = 1 };

	while (p - 1 < precision) {
		// precision is doubled every iteration
		p *= 2;
		// printf("p = %f\n", p - 1);
		p_ceil = MIN(ceil(p), precision);
		// printf("p = %f\n", p);
		// x + x * (1 - d * x)
		// d * x
		dx.size = mul(d, x, MUT_DATA(dx));
		size_t dx_point = d_point + x_point;
		/*
		printf("dx =           ");
		print_slice_bin_point(dx, dx_point);
		printf("(%zu)\n", dx_point);
		*/

		// shift 1 by d_point + x_point
		size_t shf_blocks = dx_point / BLOCK_WIDTH;
		size_t shf_bits = dx_point % BLOCK_WIDTH;
		_1_data = 1ULL << shf_bits;

		/*
		printf("%zu\n", shf);
		printf("_1 =           ");
		print_slice_bin_point(_1, shf_bits);
		printf("\n");
		*/

		// 1 - dx
		dx.size = sub_sls(_1, shf_blocks, dx, MUT_DATA(dx), &dx_sign);
		// printf("%s\n", dx_sign ? "-" : "+");

		/*
		printf("1 - dx =       ");
		print_slice_bin_point(dx, dx_point);
		printf("(%zu)\n", dx_point);
		*/

		// x * (1 - dx)
		xdx.size = mul(x, dx, MUT_DATA(xdx));
		size_t xdx_point = dx_point + x_point;
		/*
		printf("x * (1 - dx) = ");
		print_slice_bin_point(xdx, xdx_point);
		printf("(%zu)\n", xdx_point);
		*/
		if (xdx_point > p_ceil) {
			xdx.size = rshift(xdx, xdx_point - p_ceil, MUT_DATA(xdx));
			/*
			printf("Rounded      = ");
			print_slice_bin_point(xdx, p_ceil);
			printf("(%zu)\n", p_ceil);
			*/
		}
		// else printf("Unrounded\n");

		if (x_point < p_ceil) {
			x.size = lshift(x, p_ceil - x_point, MUT_DATA(x));
			x_point = p_ceil;
		}

		/*
		printf("x = ");
		print_slice_bin_point(x, p_ceil);
		printf("(%u)\n", p_ceil);
		*/
		bool x_sign;
		x.size = add_signed(x, POSITIVE, xdx, dx_sign, MUT_DATA(x), &x_sign, true);

		/*
		printf("x = ");
		print_slice_bin_point(x, p_ceil);
		printf("(%u)\n", );
		*/
	}

	/*
	print_slice(d);
	printf("\n");
	print_slice(x);
	printf("\n");
	*/

	/*
	dx.size = mul(d, x, MUT_DATA(dx));
	size_t dx_point = d_point + x_point;
	
	// shift 1 by d_point + x_point
	size_t shf_blocks = dx_point / BLOCK_WIDTH;
	size_t shf_bits = dx_point % BLOCK_WIDTH;
	_1_data = 1ULL << shf_bits;
	// print_slice_bin_point(_1, dx_point);
	// printf("\n");

	if (cmp_sls(_1, shf_blocks, dx) > 0) {
		// print_slice_bin_point(x, x_point);
		// printf("\n");
		// print_slice_bin_point(dx, dx_point);
		// printf("\n");

		x.size = add(x, ONE, MUT_DATA(x), true);
		dx.size = mul(d, x, MUT_DATA(dx));

		// print_slice_bin_point(x, x_point);
		// printf("\n");
		// print_slice_bin_point(dx, dx_point);
		// printf("\n");
		// printf("\n");

		if (cmp_sls(_1, shf_blocks, dx) > 0) {
			print_slice(d);
			printf("\n");
			print_slice(x);
			printf("\n");
			assert(0);
		}
	}
	*/

	/*
	else {
		printf("\n");
		printf("\n");
		print_slice_bin_point(dx, dx_point);
		printf("\n");
		printf("\n");
	}
	*/

	arena_pop_to(arena, restore_pos);
	return x.size;
}

static size_t div_recpr_cap(size_t num_size, size_t denum_size, size_t recpr_size, size_t precision, size_t* rem_cap) {
	assert(denum_size);
	const size_t denum_width = denum_size * BLOCK_WIDTH;
	const size_t quo_cap = num_size + recpr_size;
	const size_t shf = precision + denum_width;
	const size_t shf_quo_cap = quo_cap - shf / BLOCK_WIDTH;
	const size_t dq_cap = denum_size + shf_quo_cap;
	*rem_cap = MAX(num_size, dq_cap);
	return quo_cap;
}

// divide using reciprocal
static size_t div_recpr(Slice num, Slice denum, Slice recpr, size_t precision, Block* quo_data, Block* rem_data, size_t* rem_size) {
	assert(denum.size);
	const size_t denum_width = width(denum);
	const size_t quo_cap = num.size + recpr.size;
	const size_t shf = precision + denum_width;
	const size_t shf_quo_cap = quo_cap - shf / BLOCK_WIDTH;
	const size_t dq_cap = denum.size + shf_quo_cap;
	const size_t rem_cap = MAX(num.size, dq_cap);

	Slice quo = { .data = quo_data };
	quo.size = mul(num, recpr, quo_data);
	quo.size = rshift(quo, shf, quo_data);

	Slice dq = { .data = PUSH_ARRAY(arena, Block, dq_cap) };
	dq.size = mul(denum, quo, MUT_DATA(dq));

	Slice rem = { .data = rem_data };
	bool rem_sign;
	rem.size = sub(num, dq, rem_data, &rem_sign);
	size_t rem_width = width(rem);

	if(rem_sign == POSITIVE) {
		if (cmp(rem, denum) >= 0) {
			size_t width_diff = rem_width - denum_width;
			if (width_diff > 1) {
				WLOG_STR("BAD RECIPROCAL APPROXIMATION (remainder is %zu bits larger)", width_diff);
				quo.size = long_div(num, denum, quo_data, rem_data, &rem.size);
			} else do {
				rem.size = usub(rem, denum, rem_data);
				quo.size = add(quo, ONE, quo_data, true);
			} while (cmp(rem, denum) >= 0);
		}
	} else {
		size_t width_diff = rem_width > denum_width ? rem_width - denum_width : 0;
		if (width_diff > 1) {
			WLOG_STR("BAD RECIPROCAL APPROXIMATION (remainder is %zu bits larger and negative)", width_diff);
			quo.size = long_div(num, denum, quo_data, rem_data, rem_size);
			rem_sign = POSITIVE;
		} else do {
			bool quo_sign;
			rem.size = add_signed(rem, rem_sign, denum, POSITIVE, rem_data, &rem_sign, true);
			quo.size = sub(quo, ONE, quo_data, &quo_sign);
			assert(quo_sign == POSITIVE);
		} while(rem_sign == NEGATIVE);
		assert(cmp(rem, denum) < 0);
	}
	assert(rem_sign == POSITIVE);
	assert(cmp(rem, denum) < 0);

	*rem_size = rem.size;
	return quo.size;
}

// out = a^pow % mod
static size_t powmod(Slice a, Slice pow, Slice mod, Block* out_data) {
	assert(cmp(a, mod) < 0);
	assert(mod.size);

	Slice out = { .data = out_data, .size = 1 };
	out_data[0] = 1;

	u64 restore_pos = arena->pos;
	Slice tmp = { PUSH_ARRAY(arena, Block, 2 * mod.size) };

	for (size_t i = pow.size; i > 0; i--) {
		Block pow_block = pow.data[i - 1];
		const int bits = (i == pow.size) ?
			BLOCK_WIDTH - CLZ(pow_block) : BLOCK_WIDTH;

		for (int j = bits - 1; j >= 0; j--) {
			// tmp = out * out
			tmp.size = mul(out, out, MUT_DATA(tmp));
			// out = tmp % mod
			long_div(tmp, mod, NULL, MUT_DATA(tmp), &out.size);
			memcpy(out_data, MUT_DATA(tmp), sizeof(Block) * out.size);
			if (pow_block & ((Block)1 << j)) {
				// tmp = a * out
				tmp.size = mul(a, out, MUT_DATA(tmp));
				// out = tmp % mod
				long_div(tmp, mod, NULL, MUT_DATA(tmp), &out.size);
				memcpy(out_data, MUT_DATA(tmp), sizeof(Block) * out.size);
			}
		}
	}

	arena_pop_to(arena, restore_pos);
	return out.size;
}

/*
static size_t powmod_recpr(Slice a, Slice pow, Slice mod, Slice recpr, size_t precision, Block* out_data, Block* buffer) {
	assert(cmp(a, mod) < 0);
	assert(mod.size);
	const size_t mul_cap = mod.size * 2;
	const size_t div_cap = div_recpr_cap(2 * mod.size, mod.size, recpr.size, precision, );
	Block* const mul_buf = buffer + mul_cap;
	Block* const div_buf = buffer + mul_cap;
	Block* const quo_buf = div_buf + div_cap;
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
			div_recpr(buf, mod, recpr, precision, NULL, buffer, &out.size, div_buf);
			memcpy(out_data, buffer, sizeof(Block) * out.size);
			if (pow_block & ((Block)1 << j)) {
				// buffer = a * out
				buf.size = mul(a, out, buffer, buffer + mul_cap);
				// out = buffer % mod
				div_recpr(buf, mod, recpr, precision, NULL, buffer, &out.size, div_buf);
				memcpy(out_data, buffer, sizeof(Block) * out.size);
			}
		}
	}
}
*/

