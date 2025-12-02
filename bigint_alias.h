#pragma once
#include "bigint.h"
#include "bigint_def.h"

typedef BigInt_DataBlock  DataBlock;
typedef BigInt_CapField   CapField;
typedef BigInt_PointField PointField;
typedef BigInt_FormatSpec FormatSpec;
typedef BigInt_Slice      Slice;

#define BLOCK_WIDTH BIGINT_BLOCK_WIDTH
#define CAP_WIDTH   BIGINT_CAP_WIDTH
#define POINT_WIDTH BIGINT_POINT_WIDTH
#define MAX_CAP     BIGINT_MAX_CAP

