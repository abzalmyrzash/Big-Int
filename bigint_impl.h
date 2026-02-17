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

#define OUT(x) _Generic((x), \
	Slice: (OutSlice) { .data = (Block*)((x).data), .size = &((x).data) })

#define IN(x) _Generic((x), \
	OutSlice: (Slice) { .data = (x).data, .size = *((x).data) })

typedef struct {
	Block* buf;
	size_t bufsize;
} Context;

static thread_local mem_arena* arena;

static bool warn_bad_recpr = true;

static constexpr Block VALUE_TEN = 10;
static const Slice TEN = {
	.data = &VALUE_TEN,
	.size = 1
};
#define TEN_WIDTH 4

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
#define DDPS DEC_DIGITS_PER_SPACE

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
	#define FORMAT_DEC     "%" PRIu64

#elif BIGINT_BLOCK_WIDTH == 32
	#define ADDCARRY _addcarry_u32
	#define SUBBORROW _subborrow_u32
	typedef uint64_t BigMult;
	#define FORMAT_HEX     "%" PRIX32
	#define FORMAT_0HEX  "%08" PRIX32
	#define FORMAT_hex     "%" PRIx32
	#define FORMAT_0hex  "%08" PRIx32
	#define FORMAT_DEC     "%" PRIu32

#elif BIGINT_BLOCK_WIDTH == 16
	#define ADDCARRY _addcarry_u16
	#define SUBBORROW _subborrow_u16
	typedef uint32_t BigMult;
	#define FORMAT_HEX     "%" PRIX16
	#define FORMAT_0HEX  "%04" PRIX16
	#define FORMAT_hex     "%" PRIx16
	#define FORMAT_0hex  "%04" PRIx16
	#define FORMAT_DEC     "%" PRIu16

#elif BIGINT_BLOCK_WIDTH == 8
	#define ADDCARRY _addcarry_u8
	#define SUBBORROW _subborrow_u8
	typedef uint16_t BigMult;
	#define FORMAT_HEX     "%" PRIX8
	#define FORMAT_0HEX  "%02" PRIX8
	#define FORMAT_hex     "%" PRIx8
	#define FORMAT_0hex  "%02" PRIx8
	#define FORMAT_DEC     "%" PRIu8
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

// size in blocks
static inline size_t to_blocks(size_t width) {
	return (width + BLOCK_WIDTH - 1) / BLOCK_WIDTH;
}

static Slice split(Slice* slice, size_t size) {
	Slice second;
	if (slice->size > size) {
		second.data = slice->data + size;
		second.size = slice->size - size;
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
	const size_t n = MAX(a.size, b.size);
	const size_t m = MIN(a.size, b.size);
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
// b - first number x
// d - second number y
// n - max(x.size, y.size);
// out - filled with x * y, must be allocated outside
//       (at least (x.size + y.size) * sizeof(Block) bytes);
// @return size of out
//
// It works by splitting x into (a <<< m) + b, and y into (c <<< m) + d, where m is n/2 rounded up
// (note that b and d become different than x and y, though x and y were passed as b and d)
// then multiply [(a <<< m) + b] * [(c <<< m) + d] using formula:
// [a <<< (2 * m)] + {[(a - b)(d - c) + ac + bd] <<< m} + bd
//
// Note: (a - b)(d - c) + ac + bd is used for the middle part instead of
// the more classic (a + b)(c + d) - ac - bd
// because (a + b) and (c + d) can overflow
static size_t karatsuba(Slice b, Slice d, size_t n, Block* out) {

	if(b.size == 0 || d.size == 0) return 0;
	assert(n == MAX(b.size, d.size));

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

	// m = (n + 1) / 2, or in other words, n/2 rounded up
	const size_t m = (n + 1) / 2;
	// b is split into a and b, where a is the high part, and b is the low part
	const Slice a = split(&b, m);
	// d is split into c and d, similarly
	const Slice c = split(&d, m);
	// the lower part's size will always be >= higher part's, but roughly even

	bool ab_sign;
	bool dc_sign;
	bool abcd_sign;

	// abcd will be used to store (a - b) * (d - c),
	// then add (a * c) and (a * c) <<< m
	// then add (b * d)
	// shift the whole thing by m blocks,
	// and finally add (b * d) again
	Slice abcd      = { .data = out };
	
	size_t restore_pos = arena->pos;
	// max cap for buffer is 2 * m blocks, as evidenced below
	Block* const buffer = PUSH_ARRAY(arena, Block, 2 * m);

	// first we use buffer to store (a - b) and (d - c)
	Slice a_minus_b = { .data = buffer };     // max cap = m, since both a and b's sizes are <= m
	Slice d_minus_c = { .data = buffer + m }; // max cap = m, since both d and c's sizes are <= m
	// (a - b) and (d - c) together use 2 * m

	// after multiplying them (storing the result in abcd),
	// the buffer will be free to be used by (a * c)
	Slice a_times_c = { .data = buffer };     // max cap = 2 * (n - m), since a.size and c.size are <= (n - m)
											  // and multiplying a and c requires a.size + c.size
											  // NOTE: (n - m) <= m, since m = (n + 1) / 2

	// after we're done using (a * c),
	// the buffer will be free to be used by (b * d)
	Slice b_times_d = { .data = buffer };     // max cap = 2 * m, since b.size and d.size are <= m

	// as we can see, we only ever use 2 * m blocks for buffer (for one call of this function)

	a_minus_b.size = sub(a, b, MUT_DATA(a_minus_b), &ab_sign);
	assert(a_minus_b.size <= m);
	assert(a_minus_b.size == size_without_zeros(a_minus_b.data, a_minus_b.size));

	d_minus_c.size = sub(d, c, MUT_DATA(d_minus_c), &dc_sign);
	assert(d_minus_c.size <= m);
	assert(d_minus_c.size == size_without_zeros(d_minus_c.data, d_minus_c.size));

	abcd_sign = ab_sign ^ dc_sign;

	// abcd = (a-b)(d-c)
	abcd.size = mul(a_minus_b, d_minus_c, MUT_DATA(abcd));
	assert(abcd.size <= 2 * m);
	assert(abcd.size == size_without_zeros(abcd.data, abcd.size));

	a_times_c.size = mul(a, c, MUT_DATA(a_times_c));
	assert(a_times_c.size <= 2 * (n - m));
	assert(a_times_c.size == size_without_zeros(a_times_c.data, a_times_c.size));

	// abcd = abcd + ac
	abcd.size = add_signed(abcd, abcd_sign, a_times_c, POSITIVE, MUT_DATA(abcd), &abcd_sign, true);
	assert(abcd.size == size_without_zeros(abcd.data, abcd.size));
	// abcd = (a_times_c <<< m) + abcd
	abcd.size = add_sls_signed(a_times_c, m, POSITIVE, abcd, abcd_sign, MUT_DATA(abcd), &abcd_sign, true);
	assert(abcd.size == size_without_zeros(abcd.data, abcd.size));

	b_times_d.size = mul(b, d, MUT_DATA(b_times_d));
	assert(b_times_d.size <= 2 * m);
	assert(b_times_d.size == size_without_zeros(b_times_d.data, b_times_d.size));

	// abcd = abcd + bd
	abcd.size = add_signed(abcd, abcd_sign, b_times_d, POSITIVE, MUT_DATA(abcd), &abcd_sign, true);
	assert(abcd.size == size_without_zeros(abcd.data, abcd.size));
	// move abcd from out to out + m, basically leftshifting abcd <<< m
	memmove(out + m, out, abcd.size * sizeof(Block));
	abcd.data += m;
	// (abcd <<< m) + bd
	size_t size = add_sls_signed(abcd, m, abcd_sign, b_times_d, 0, out, &abcd_sign, true);
	arena_pop_to(arena, restore_pos);
	assert(size == size_without_zeros(out, size));
	return size;
}

static size_t move(Slice z, Block* out) {
	memmove(out, z.data, z.size * sizeof(Block));
	return z.size;
}

static size_t copy(Slice z, Block* out) {
	memcpy(out, z.data, z.size * sizeof(Block));
	return z.size;
}

// super-left-shifts z, i.e. left-shifts by number of blocks rather than bits
// this operation is denoted as <<< (as opposed to << for left-shift by bits)
static size_t super_lshift(Slice z, size_t shift, Block* out) {
	if (z.size == 0) return 0;
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

	while (i > 0) {
		i--;
		hi = z.data[i] & HIGH_BITS;
		lo = z.data[i] & LOW_BITS;
		out[i + shift_blocks] = lo << shift_bits;
		out[i + shift_blocks + 1] |= hi >> (BLOCK_WIDTH - shift_bits);
	}

	memset(out, 0, shift_blocks * sizeof(Block));
	return out_size;
}

static size_t lshift_1(Slice z, Block* out) {
	if (z.size == 0) return 0;

	const Block LOW_BITS = (1ULL << (BLOCK_WIDTH - 1)) - 1;
	const Block HIGH_BITS = ~0ULL - LOW_BITS;

	size_t out_size = z.size;

	size_t i = z.size - 1;
	Block hi, lo;

	hi = z.data[i] & HIGH_BITS;
	lo = z.data[i] & LOW_BITS;
	if (hi) {
		out_size++;
		out[i + 1] = hi >> (BLOCK_WIDTH - 1);
	}
	out[i] = lo << 1;

	if (i--) while (1) {
		hi = z.data[i] & HIGH_BITS;
		lo = z.data[i] & LOW_BITS;
		out[i] = lo << 1;
		out[i + 1] |= hi >> (BLOCK_WIDTH - 1);
		if (i == 0) break;
		i--;
	}

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

static inline bool bit_get(const Block* z, size_t bitpos) {
	return z[bitpos / BLOCK_WIDTH] & (1ULL << (bitpos % BLOCK_WIDTH));
}

static inline void bit_set(Block* z, size_t bitpos) {
	z[bitpos / BLOCK_WIDTH] |= (1ULL << (bitpos % BLOCK_WIDTH));
}

static inline void bit_unset(Block* z, size_t bitpos) {
	z[bitpos / BLOCK_WIDTH] &= ~(1ULL << (bitpos % BLOCK_WIDTH));
}

static inline void bit_set_to(Block* z, size_t bitpos, bool val) {
	z[bitpos / BLOCK_WIDTH] |= ((uint64_t)val << (bitpos % BLOCK_WIDTH));
	z[bitpos / BLOCK_WIDTH] &= ~((uint64_t)!val << (bitpos % BLOCK_WIDTH));
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
		*rem_size = move(a, rem_data);
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
		*rem_size = move(a, rem_data);
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
		.size = move(a, rem_data)
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
			c.size = move(b, MUT_DATA(c));
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

	bitpos = (i64)size_diff * BLOCK_WIDTH - (i64)clz_diff;

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

	assert(quo_size == size_without_zeros(quo, quo_size));
	assert(rem.size == size_without_zeros(rem.data, rem.size));
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

	size_t restore_pos = arena->pos;

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

	if (x_point > p_ceil) {
		x.size = rshift(x, x_point - p_ceil, MUT_DATA(x));
		x_point = p_ceil;
	}

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
		p_ceil = MIN(ceil(p), precision);

		// x + x * (1 - d * x)
		// d * x
		dx.size = mul(d, x, MUT_DATA(dx));
		size_t dx_point = d_point + x_point;

		// shift 1 by d_point + x_point
		size_t shf_blocks = dx_point / BLOCK_WIDTH;
		size_t shf_bits = dx_point % BLOCK_WIDTH;
		_1_data = 1ULL << shf_bits;

		// 1 - dx
		dx.size = sub_sls(_1, shf_blocks, dx, MUT_DATA(dx), &dx_sign);

		// x * (1 - dx)
		xdx.size = mul(x, dx, MUT_DATA(xdx));
		size_t xdx_point = dx_point + x_point;

		// shift to p_ceil precision
		if (xdx_point > p_ceil) {
			xdx.size = rshift(xdx, xdx_point - p_ceil, MUT_DATA(xdx));
		}

		if (x_point < p_ceil) {
			x.size = lshift(x, p_ceil - x_point, MUT_DATA(x));
			x_point = p_ceil;
		}

		// x + x * (1 - dx)
		bool x_sign;
		x.size = add_signed(x, POSITIVE, xdx, dx_sign, MUT_DATA(x), &x_sign, true);
	}

	arena_pop_to(arena, restore_pos);
	assert(x.size == size_without_zeros(x.data, x.size));
	return x.size;
}

static size_t div_recpr_cap(size_t num_size, size_t denum_size, size_t recpr_size, size_t precision, size_t* rem_cap) {
	assert(denum_size);
	const size_t denum_width = denum_size * BLOCK_WIDTH;
	const size_t quo_cap = num_size + recpr_size;
	const size_t shf = (precision + denum_width);
	const size_t shf_quo_cap = (quo_cap < shf / BLOCK_WIDTH) ? 0 : (quo_cap - shf / BLOCK_WIDTH + 1);
	const size_t dq_cap = denum_size + shf_quo_cap;
	*rem_cap = MAX(num_size, dq_cap);
	return quo_cap;
}

static int cnt_div = 0;
static int cnt_bad_recpr = 0;

#define RECIPROCAL_TOLERANCE 4

// divide using reciprocal
static size_t div_recpr(Slice num, Slice denum, Slice recpr, size_t precision, Block* quo_data, Block* rem_data, size_t* rem_size) {
	assert(denum.size);
	const size_t denum_width = width(denum);
	const size_t quo_cap = num.size + recpr.size;
	const size_t shf = (precision + denum_width);
	const size_t shf_quo_cap = (quo_cap < shf / BLOCK_WIDTH) ? 0 : (quo_cap - shf / BLOCK_WIDTH + 1);
	const size_t dq_cap = denum.size + shf_quo_cap;
	const size_t rem_cap = MAX(num.size, dq_cap);

	size_t restore_pos = arena->pos;

	Slice quo = { .data = quo_data };
	quo.size = mul(num, recpr, quo_data);
	quo.size = rshift(quo, shf, quo_data);

	Slice dq = { .data = PUSH_ARRAY(arena, Block, dq_cap) };
	dq.size = mul(denum, quo, MUT_DATA(dq));

	Slice rem = { .data = rem_data };
	bool rem_sign;
	rem.size = sub(num, dq, rem_data, &rem_sign);
	size_t rem_width = width(rem);
	cnt_div++;

	if(rem_sign == POSITIVE) {
		if (cmp(rem, denum) >= 0) {
			size_t width_diff = rem_width - denum_width;
			if (width_diff > RECIPROCAL_TOLERANCE) {
				cnt_bad_recpr++;
				if (warn_bad_recpr) WLOG_STR("BAD RECIPROCAL APPROXIMATION (remainder is %zu bits larger)", width_diff);
				quo.size = long_div(num, denum, quo_data, rem_data, &rem.size);
			} else do {
				rem.size = usub(rem, denum, rem_data);
				quo.size = add(quo, ONE, quo_data, true);
			} while (cmp(rem, denum) >= 0);
		}
	} else {
		size_t width_diff = rem_width > denum_width ? rem_width - denum_width : 0;
		if (width_diff > RECIPROCAL_TOLERANCE) {
			cnt_bad_recpr++;
			if (warn_bad_recpr) WLOG_STR("BAD RECIPROCAL APPROXIMATION (remainder is %zu bits larger and negative)", width_diff);
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

	arena_pop_to(arena, restore_pos);

	*rem_size = rem.size;
	assert(quo.size == size_without_zeros(quo.data, quo.size));
	assert(rem.size == size_without_zeros(rem.data, rem.size));
	return quo.size;
}

static size_t power(Slice a, size_t exp, Block* out_data) {
	if (exp == 0) {
		out_data[0] = 1;
		return 1;
	}

	Slice out = { .data = out_data, .size = 1 };
	out_data[0] = 1;

	size_t restore_pos = arena->pos;
	Block* tmp = PUSH_ARRAY(arena, Block, to_blocks(exp * width(a)) + 1);

	const int bits = BLOCK_WIDTH - CLZ(exp);

	for (int j = bits - 1; j >= 0; j--) {
		// out = out * out
		out.size = mul(out, out, tmp);
		memcpy(out_data, tmp, out.size * sizeof(Block));
		if (exp & ((Block)1 << j)) {
			// out = a * out
			out.size = mul(a, out, tmp);
			memcpy(out_data, tmp, out.size * sizeof(Block));
		}
	}

	arena_pop_to(arena, restore_pos);
	assert(out.size == size_without_zeros(out.data, out.size));
	return out.size;
}

// out = a^exp % mod
static size_t powmod(Slice a, Slice exp, Slice mod, Block* out_data) {
	assert(mod.size);
	assert(out_data != a.data);

	int cmp_a_mod = cmp(a, mod);
	if (a.size == 0 || cmp_a_mod == 0) {
		if (exp.size > 0) return 0;
		out_data[0] = 1;
		return 1;
	}

	size_t restore_pos = arena->pos;

	if (cmp_a_mod > 0) {
		Block* rem = PUSH_ARRAY(arena, Block, a.size);
		long_div(a, mod, NULL, rem, &a.size);
		if (a.size == 0) {
			arena_pop_to(arena, restore_pos);
			if (exp.size > 0) return 0;
			out_data[0] = 1;
			return 1;
		}
		a.data = rem;
	}

	Slice out = { .data = out_data, .size = 1 };
	out_data[0] = 1;

	Slice tmp = { PUSH_ARRAY(arena, Block, 2 * mod.size) };

	for (size_t i = exp.size; i > 0; i--) {
		Block exp_block = exp.data[i - 1];
		const int bits = (i == exp.size) ?
			BLOCK_WIDTH - CLZ(exp_block) : BLOCK_WIDTH;

		for (int j = bits - 1; j >= 0; j--) {
			// tmp = out * out
			tmp.size = mul(out, out, MUT_DATA(tmp));
			// tmp = tmp % mod
			// we don't pass out_data to div directly, as rem_data could require more space than out
			long_div(tmp, mod, NULL, MUT_DATA(tmp), &out.size);
			// copy tmp to out 
			memcpy(out_data, tmp.data, sizeof(Block) * out.size);
			if (exp_block & ((Block)1 << j)) {
				// tmp = a * out
				tmp.size = mul(a, out, MUT_DATA(tmp));
				// tmp = tmp % mod, copy tmp to out
				long_div(tmp, mod, NULL, MUT_DATA(tmp), &out.size);
				memcpy(out_data, tmp.data, sizeof(Block) * out.size);
			}
		}
	}

	arena_pop_to(arena, restore_pos);
	assert(out.size == size_without_zeros(out.data, out.size));
	return out.size;
}

static size_t powmod_recpr(Slice a, Slice exp, Slice mod, Block* out_data) {
	assert(mod.size);
	assert(out_data != a.data);

	int cmp_a_mod = cmp(a, mod);
	if (a.size == 0 || cmp_a_mod == 0) {
		if (exp.size > 0) return 0;
		out_data[0] = 1;
		return 1;
	}

	size_t restore_pos = arena->pos;

	size_t mod_width = width(mod);

	if (cmp_a_mod > 0) {
		size_t a_width = width(a);
		size_t precision = MAX(a_width, mod_width);
		size_t recpr_cap = newton_reciprocal_cap(precision);
		Slice recpr = { PUSH_ARRAY(arena, Block, recpr_cap) };
		recpr.size = newton_reciprocal(mod, precision, MUT_DATA(recpr));

		size_t quo_cap, rem_cap;
		quo_cap = div_recpr_cap(a.size, mod.size, recpr.size, precision, out_data);
		Block* rem = PUSH_ARRAY(arena, Block, rem_cap);
		Block* quo = PUSH_ARRAY(arena, Block, quo_cap);
		div_recpr(a, mod, recpr, precision, quo, rem, &a.size);
		if (a.size == 0) {
			arena_pop_to(arena, restore_pos);
			if (exp.size > 0) return 0;
			out_data[0] = 1;
			return 1;
		}
		POP_ARRAY(arena, Block, quo_cap);
		a.data = rem;
	}

	size_t precision = mod_width * 2;
	size_t recpr_cap = newton_reciprocal_cap(precision);
	Slice recpr = { PUSH_ARRAY(arena, Block, recpr_cap) };
	recpr.size = newton_reciprocal(mod, precision, MUT_DATA(recpr));

	Slice out = { .data = out_data, .size = 1 };
	out_data[0] = 1;

	size_t num_cap = mod.size * 2;
	size_t quo_cap;
	size_t rem_cap;

	quo_cap = div_recpr_cap(num_cap, mod.size, recpr.size, precision, &rem_cap);

	Slice num = { .data = PUSH_ARRAY(arena, Block, num_cap) };
	Block* quo = PUSH_ARRAY(arena, Block, quo_cap);
	Block* rem = PUSH_ARRAY(arena, Block, rem_cap);

	for (size_t i = exp.size; i > 0; i--) {
		Block pow_block = exp.data[i - 1];
		const int bits = (i == exp.size) ?
			BLOCK_WIDTH - CLZ(pow_block) : BLOCK_WIDTH;
		for (int j = bits - 1; j >= 0; j--) {
			// num = out * out
			num.size = mul(out, out, MUT_DATA(num));
			// rem = num % mod
			// we don't pass out_data to div directly, as rem_data could require more space than out
			div_recpr(num, mod, recpr, precision, quo, rem, &out.size);
			// copy rem to out
			memcpy(out_data, rem, sizeof(Block) * out.size);
			if (pow_block & ((Block)1 << j)) {
				// num = a * out
				num.size = mul(a, out, MUT_DATA(num));
				// rem = num % mod, copy rem to out
				div_recpr(num, mod, recpr, precision, quo, rem, &out.size);
				memcpy(out_data, rem, sizeof(Block) * out.size);
			}
		}
	}
	arena_pop_to(arena, restore_pos);
	assert(out.size == size_without_zeros(out.data, out.size));
	return out.size;
}

/*
0000 0000 1011
0000 0001 0110
0000 0010 1100
0000 0101 1000
0000 1000 1000
0001 0001 0000
*/

#define GET_BIT(x, bit) (x)[(bit) / BLOCK_WIDTH] & (1ULL << ((bit) % BLOCK_WIDTH))
#define SET_BIT(x, bit) ((Block*)(x))[(bit) / BLOCK_WIDTH] |= (1ULL << ((bit) % BLOCK_WIDTH))
#define ASSIGN_BIT(x, bit, val) ((Block*)(x))[(bit) / BLOCK_WIDTH] |= ((Block)(bool)(val) << ((bit) % BLOCK_WIDTH))

static size_t decimal_width(size_t bit_width) {
	return ceil(bit_width * (LOG10_2 + 1e-6));
}

static size_t divmod10(Slice x, char* str) {
	if (x.size == 0) {
		return 0;
	}

	const size_t q_cap = x.size;
	const size_t tmp_cap = x.size;

	size_t restore_pos = arena->pos;
	Block* const q_data = PUSH_ARRAY(arena, Block, q_cap);
	Block* const tmp_data = PUSH_ARRAY(arena, Block, tmp_cap);
	Slice tmp = {
		.data = tmp_data,
		.size = move(x, tmp_data)
	};

	USmallInt rem;
	size_t rem_size;

	size_t cnt = 0;
	while (tmp.size > 0) {
		tmp.size = long_div(tmp, TEN, q_data, tmp_data, &rem_size);
		assert((rem_size == 1 && tmp_data[0] < 10) || rem_size == 0);
		rem = rem_size ? tmp_data[0] : 0;
		str[cnt++] = rem + '0';
		memcpy(tmp_data, q_data, tmp.size * sizeof(Block));
	}

	for (size_t i = 0; i < cnt / 2; i++) {
		char tmp = str[i];
		str[i] = str[cnt - i - 1];
		str[cnt - i - 1] = tmp;
	}
	str[cnt] = '\0';

	arena_pop_to(arena, restore_pos);
	return cnt;
}

static size_t divmod10_recpr(Slice x, char* str) {
	if (x.size == 0) {
		return 0;
	}
	size_t restore_pos = arena->pos;

	size_t precision = MAX(width(x), TEN_WIDTH);
	size_t new_p;
	size_t recpr_cap = newton_reciprocal_cap(precision);

	Slice recpr = { .data = PUSH_ARRAY(arena, Block, recpr_cap) };
	recpr.size = newton_reciprocal(TEN, precision, MUT_DATA(recpr));

	size_t quo_cap = x.size;
	size_t num_cap = x.size;
	size_t rem_cap;

	quo_cap = div_recpr_cap(x.size, TEN.size, recpr.size, precision, &rem_cap);

	Block* const quo_data = PUSH_ARRAY(arena, Block, quo_cap);
	Block* const num_data = PUSH_ARRAY(arena, Block, num_cap);
	Slice num = {
		.data = num_data,
		.size = move(x, num_data)
	};

	Block* rem_data = PUSH_ARRAY(arena, Block, rem_cap);
	size_t rem_size;

	size_t cnt = 0;
	while (num.size > 0) {
		num.size = div_recpr(num, TEN, recpr, precision, quo_data, rem_data, &rem_size);
		assert((rem_size == 1 && rem_data[0] < 10) || rem_size == 0);
		Block rem = rem_size ? rem_data[0] : 0;
		str[cnt++] = rem + '0';
		memcpy(num_data, quo_data, num.size * sizeof(Block));
		new_p = MAX(width(num), TEN_WIDTH);
		assert(new_p <= precision);
		recpr.size = rshift(recpr, precision - new_p, MUT_DATA(recpr));
		precision = new_p;
	}

	for (size_t i = 0; i < cnt / 2; i++) {
		char tmp = str[i];
		str[i] = str[cnt - i - 1];
		str[cnt - i - 1] = tmp;
	}
	str[cnt] = '\0';

	arena_pop_to(arena, restore_pos);
	return cnt;
}

static size_t double_dabble(Slice x, char* str) {
	if (x.size == 0) {
		return 0;
	}
	size_t restore_pos = arena->pos;

	const size_t n = width(x);
	const size_t buf_max_width = 4 * decimal_width(n);
	const size_t buf_max_size = ceil_div(buf_max_width, BLOCK_WIDTH);
	Slice buf = { .data = PUSH_ARRAY_ZERO(arena, Block, buf_max_size) };
	MUT(buf.data)[0] = 1;
	size_t width = 1;
	buf.size = 1;
	size_t bcd_cnt = 1;
	constexpr size_t BCD_PER_BLOCK = BLOCK_WIDTH / 4;
	Block mask = 15;

	for (size_t i = 2; i <= n; i++) {
		for (size_t j = 0; j < buf.size; j++) {
			int shift = 0;
			for (int k = 0; k < BCD_PER_BLOCK; k++) {
				Block bcd = (buf.data[j] & (mask << shift)) >> shift;
				if (bcd >= 5) {
					MUT(buf.data)[j] += 3ULL << shift;
				}
				shift += 4;
			}
		}
		size_t old_size = buf.size;
		buf.size = lshift_1(buf, MUT_DATA(buf));
		if (CLZ(buf.data[buf.size - 1]) % 4 == 3) {
			bcd_cnt++;
		}
		ASSIGN_BIT(buf.data, 0, GET_BIT(x.data, n - i));
	}

	int shift = 0;
	for (size_t i = 0; i < bcd_cnt; i++) {
		str[bcd_cnt - i - 1] = '0' + ((buf.data[i / BCD_PER_BLOCK] & (mask << shift)) >> shift);
		shift += 4;
		shift %= BLOCK_WIDTH;
	}

	arena_pop_to(arena, restore_pos);
	return bcd_cnt;
}

size_t decimal_split(Slice x, char* str) {
	if (x.size == 0) {
		return 0;
	}
	if (x.size == 1) {
		static char tmp[MAX_DEC_INT_DIGITS(Block)];
		int len = sprintf(tmp, FORMAT_DEC, x.data[0]);
		memcpy(str, tmp, sizeof(char) * len);
		return len;
	}
	size_t restore_pos = arena->pos;

	size_t x_width = width(x);
	size_t dec_width = decimal_width(x_width);
	size_t half_dw = dec_width / 2;
	Slice ten_pow_hdw = { PUSH_ARRAY(arena, Block, to_blocks(half_dw * TEN_WIDTH)) };
	ten_pow_hdw.size = power(TEN, half_dw, MUT_DATA(ten_pow_hdw));
	size_t precision = MAX(x_width, width(ten_pow_hdw));
	Slice recpr = { PUSH_ARRAY(arena, Block, newton_reciprocal_cap(precision)) };
	recpr.size = newton_reciprocal(ten_pow_hdw, precision, MUT_DATA(recpr));

	size_t quo_cap, rem_cap;
	quo_cap = div_recpr_cap(x.size, ten_pow_hdw.size, recpr.size, precision, &rem_cap);
	Slice quo = { PUSH_ARRAY(arena, Block, quo_cap) };
	Slice rem = { PUSH_ARRAY(arena, Block, rem_cap) };
	quo.size = div_recpr(x, ten_pow_hdw, recpr, precision, MUT_DATA(quo), MUT_DATA(rem), &rem.size);
	size_t quo_dw = decimal_split(quo, str);
	size_t rem_dw = decimal_split(rem, str + quo_dw);

	arena_pop_to(arena, restore_pos);

	if (quo_dw > 0) {
		if (rem_dw < half_dw) {
			memmove(str + quo_dw + half_dw - rem_dw, str + quo_dw, sizeof(char) * rem_dw);
			memset(str + quo_dw, '0', sizeof(char) * (half_dw - rem_dw));
		}
		return quo_dw + half_dw;
	} else return rem_dw;
}

//            1          2          3           4
// 01234567 89012345 67890123 45678901 23456789 01234567
//                ^                              ^
// bitpos = 14, width = 28, block_width = 16
// byte_offset = bitpos / 8 = 1, bit_offset = bitpos % 8 = 6
// len = width / 8 = 3, last_bits = width % 8 = 4
// block = 0, shift = 0
// for i = 0 to len {
//     byte = (high (8 - bit_offset = 2) bits from [i + byte_offset] >> bit_offset)
//        | (low (bit_offset = 6) bits [i + byte_offset + 1] << (8 - bit_offset))
//     dst[block] = byte << shift
//     shift += 8
//     if shift reaches block_width, increment block, set shift = 0
// }
// 45678901 23456789 | 01234567
// if (last_bits > 0) {
//     byte = high 2 bits from [len + byte_offset]
//     if (rem_bits > (8 - bit_offset))
//         byte |= low 6 bits from [len + byte_offset + 1] << 2
//     byte &= mask for lowest rem_bits
//     dst[block] = byte << shift
//     block++
// }
// 45678901 23456789 | 01234567 8901


static size_t from_little_endian(const u8* src_bytes, size_t bitpos, size_t width, Block* dst_blocks) {
	size_t max_size = ceil_div(width, BLOCK_WIDTH);
	memset(dst_blocks, 0, sizeof(Block) * max_size);

	const size_t byte_offset = bitpos / CHAR_BIT;
	const size_t bit_offset  = bitpos % CHAR_BIT;
	const size_t bit_remwid  = CHAR_BIT - bit_offset;
	const size_t len         = width / CHAR_BIT;
	const size_t last_bits   = width % CHAR_BIT;

	const u8 LOW_BITS = (u8)(1U << bit_offset) - 1;
	const u8 HIGH_BITS = (u8)~0U - LOW_BITS;

	u8 byte;
	size_t block = 0;
	size_t shift = 0;
	for (size_t i = 0; i < len; i++) {
		byte = (src_bytes[i + byte_offset] & HIGH_BITS) >> bit_offset;
		if (LOW_BITS) byte |= (src_bytes[i + byte_offset + 1] & LOW_BITS) << bit_remwid;
		assert(block < max_size);
		dst_blocks[block] |= (Block)byte << shift;
		shift += CHAR_BIT;
		if (shift == BLOCK_WIDTH) {
			block++;
			shift = 0;
		}
	}

	if (last_bits > 0) {
		byte = (src_bytes[len + byte_offset] & HIGH_BITS) >> bit_offset;
		if (last_bits > bit_remwid) {
			byte |= (src_bytes[len + byte_offset + 1] & LOW_BITS) << bit_remwid;
		}
		byte &= (1 << last_bits) - 1;
		assert(block < max_size);
		dst_blocks[block] |= (Block)byte << shift;
	}
	fflush(stdout);

	return size_without_zeros(dst_blocks, max_size);
}

// 45678901 23456789 | 01234567 8901
//            1          2          3           4
// 01234567 89012345 67890123 45678901 23456789 01234567
//                ^                              ^
// bitpos = 14, width = 28, block_width = 16
// byte_offset = bitpos / 8 = 1, bit_offset = bitpos % 8 = 6
// len = width / 8 = 3, rem_bits = width % 8 = 4
// [byte_offset] &= low_mask (8th-13th bits)
// 01234567 890123xx
// block = 0, shift = 0
// byte_mask = 255
// for i = 0 to len:
//     byte = [block] & (byte_mask << shift) >> shift
//     shift += 8
//     dst[byte_offset + i] |= (byte & low_mask) << bit_offset
//     dst[byte_offset + i + 1] = (byte & high_mask) >> (8 - bit_offset)
//     if (shift == block_width) {
//          shift = 0
//          block++
//     }
//  01234567 89012345 67890123 45678901 234567xx xx
//  add 8901
//  if (last_bits > 0) {
//      byte = [block] & (byte_mask << shift) >> shift
//      dst[byte_offset + len] |= (byte & low_mask) << bit_offset
//      dst[byte_offset + len + 1] &= ~0 - (1 << bit_offset);
//  }

static void to_little_endian(const Block* src_blocks, size_t width, u8* dst_bytes, size_t bitpos) {
	const size_t byte_offset = bitpos / CHAR_BIT;
	const size_t bit_offset  = bitpos % CHAR_BIT;
	const size_t bit_remwid  = CHAR_BIT - bit_offset;
	const size_t len         = width / CHAR_BIT;
	const size_t last_bits   = width % CHAR_BIT;

	const u8 LOW_BITS  = (u8)((1 << bit_remwid) - 1);
	const u8 HIGH_BITS = (u8)~0 - LOW_BITS;

	u8 byte;
	size_t block = 0;
	int shift = 0;
	constexpr Block BYTE_MASK = 255;
	for (size_t i = 0; i < len; i++) {
		byte = (src_blocks[block] & (BYTE_MASK << shift)) >> shift;
		dst_bytes[byte_offset + i] |= (byte & LOW_BITS) << bit_offset;
		if (HIGH_BITS) dst_bytes[byte_offset + i + 1] = (byte & HIGH_BITS) >> bit_remwid;
		shift += 8;
		if (shift == BLOCK_WIDTH) {
			shift = 0;
			block++;
		}
	}

	if (last_bits > 0) {
		byte = (src_blocks[block] & (BYTE_MASK << shift)) >> shift;
		dst_bytes[byte_offset + len] |= (byte & LOW_BITS) << bit_offset;
		if (last_bits > bit_remwid) {
			dst_bytes[byte_offset + len + 1] = (byte & HIGH_BITS) >> bit_remwid;
		}
	}
}

