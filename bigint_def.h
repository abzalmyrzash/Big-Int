#pragma once
#include "bigint_params.h"

// size:  length of number in datablocks
// cap:   number of datablocks allocated
// point: position of the fractional point
// sign:  0 means positive, 1 means negative
struct bigint {
	BigInt_CapField   size;
	BigInt_CapField   cap;
#ifdef BIGINT_SPLIT_POINT_AND_SIGN
	BigInt_PointField point : BIGINT_POINT_WIDTH - 1;
	BigInt_PointField sign  : 1;
#else
	BigInt_PointField point;
	bool sign;
#endif
	BigInt_DataBlock  data[];
}
#ifdef BIGINT_PACKED
	__attribute__((packed))
#endif
;

typedef struct {
	const BigInt_DataBlock* data;
	BigInt_CapField size;
} BigInt_Slice;

