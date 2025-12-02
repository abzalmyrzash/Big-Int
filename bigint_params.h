#pragma once
#include <stdint.h>

// width (size in bits) of fields in struct bigint
#define BIGINT_BLOCK_WIDTH 64  // width of one data block
#define BIGINT_CAP_WIDTH   32  // width of capacity field; also equal to size width
#define BIGINT_POINT_WIDTH 32  // width of point field

// if defined one bit of point field will be used for sign (like a bitfield)
// otherwise sign will be separate bool field
#define BIGINT_SPLIT_POINT_AND_SIGN 

// maximum capacity
#define BIGINT_MAX_CAP     (uint64_t)~0ULL >> (64 - BIGINT_CAP_WIDTH)

#define BIGINT_KARATSUBA_THRESHOLD 9

// pack the bigint struct
// #define BIGINT_PACKED 

#if BIGINT_BLOCK_WIDTH == 64
	typedef  int64_t SmallInt;
	typedef uint64_t USmallInt;
	typedef uint64_t BigInt_DataBlock;
#elif BIGINT_BLOCK_WIDTH == 32
	typedef  int32_t SmallInt;
	typedef uint32_t USmallInt;
	typedef uint32_t BigInt_DataBlock;
#elif BIGINT_BLOCK_WIDTH == 16
	typedef  int16_t SmallInt;
	typedef uint16_t USmallInt;
	typedef uint16_t BigInt_DataBlock;
#elif BIGINT_BLOCK_WIDTH == 8
	typedef  int8_t  SmallInt;
	typedef uint8_t  USmallInt;
	typedef uint8_t  BigInt_DataBlock;
#else
	#error "Invalid BIGINT_BLOCK_WIDTH (only 8, 16, 32, 64 are allowed)"
#endif

#if BIGINT_CAP_WIDTH == 64
	typedef uint64_t BigInt_CapField;
#elif BIGINT_CAP_WIDTH == 32
	typedef uint32_t BigInt_CapField;
#elif BIGINT_CAP_WIDTH == 16
	typedef uint16_t BigInt_CapField;
#elif BIGINT_CAP_WIDTH == 8
	typedef uint8_t  BigInt_CapField;
#else
	#error "Invalid BIGINT_BLOCK_WIDTH (only 8, 16, 32, 64 are allowed)"
#endif

#if BIGINT_POINT_WIDTH == 64
	typedef uint64_t BigInt_PointField;
#elif BIGINT_POINT_WIDTH == 32
	typedef uint32_t BigInt_PointField;
#elif BIGINT_POINT_WIDTH == 16
	typedef uint16_t BigInt_PointField;
#elif BIGINT_POINT_WIDTH == 8
	typedef uint8_t  BigInt_PointField;
#else
	#error "Invalid BIGINT_BLOCK_WIDTH (only 8, 16, 32, 64 are allowed)"
#endif

