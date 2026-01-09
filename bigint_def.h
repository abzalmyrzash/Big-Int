#pragma once
#include "bigint_params.h"

// size:  length of number in datablocks
// cap:   number of datablocks allocated
// sign:  0 means positive, 1 means negative
struct bigint {
	BigInt_CapField   cap;
#ifdef BIGINT_SPLIT_SIZE_AND_SIGN
	BigInt_CapField   size  : BIGINT_CAP_WIDTH - 1;
	BigInt_CapField   sign  : 1;
#else
	BigInt_CapField   size;
	bool sign;
#endif
	BigInt_Block  data[];
}
#ifdef BIGINT_PACKED
	__attribute__((packed))
#endif
;

typedef struct {
	const BigInt_Block* data;
	size_t size;
} BigInt_Slice;

