#include "bigint_rand.h"
#include "bigint_alias.h"
#include "bigint_impl_basic.h"
#include <stdlib.h>
#include "utils.h"

static constexpr int RAND_WIDTH = 32 - __builtin_clzg((uint32_t)RAND_MAX);

BigInt bigint_rand(BigInt* z_ptr, size_t bits) {
	const size_t blocks = bits / BLOCK_WIDTH;
	const int final_bits = bits % BLOCK_WIDTH;
	const size_t cap = blocks + (final_bits > 0);
	BigInt z = bigint_reserve(z_ptr, cap, CLEAR_DATA);
	int rand_num = 0;
	int rem_bits = 0;
	for (size_t i = 0; i < blocks; i++) {
		z->data[i] = rand_num;
		int j;
		for (j = rem_bits; j < BLOCK_WIDTH; j += RAND_WIDTH) {
			rand_num = rand();
			// rand_num &= MIN(RAND_WIDTH, BLOCK_WIDTH - j);
			z->data[i] |= (Block)rand_num << j;
		}
		rem_bits = j - BLOCK_WIDTH;
		rand_num >>= RAND_WIDTH - rem_bits;
	}
	if (final_bits > 0) {
		const Block final_bits_mask = ((Block)1 << final_bits) - 1;
		z->data[blocks] = rand_num;
		for (int i = rem_bits; i < final_bits; i += RAND_WIDTH) {
			rand_num = rand();
			// rand_num &= MIN(RAND_WIDTH, BLOCK_WIDTH - i);
			z->data[blocks] |= (Block)rand_num << i;
		}
		z->data[blocks] &= final_bits_mask;
	}
	z->size = size_without_zeros(z->data, cap);
	return z;
}

