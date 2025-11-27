#pragma once
#include "bigint_params.h"
#include <stdio.h>
#include <stdbool.h>

extern int bigint_errno;
#define BIGINT_ERR_ALLOC_FAIL  -1
#define BIGINT_ERR_SUB_NEG     -2
#define BIGINT_ERR_DIV_BY_ZERO -3

typedef struct bigint_context* BigInt_Context;

typedef struct bigint* BigInt;
typedef const struct bigint* ConstBigInt;
typedef struct {
	uint8_t base;
	bool is_unsigned;
	bool add_plus_sign;
	bool add_prefix;
	bool leading_zeros;
	bool add_spaces;
	bool uppercase;
} BigInt_FormatSpec;

void bigint_init();
void bigint_finish();

BigInt bigint_alloc(size_t cap);
BigInt bigint_zalloc(size_t cap);
BigInt bigint_realloc(BigInt* z, size_t cap, bool keep_data);

void bigint_free(BigInt z);
void bigint_structinfo();
void bigint_memstat();

size_t bigint_size(ConstBigInt z);
size_t bigint_cap(ConstBigInt z);
size_t bigint_point(ConstBigInt z);
bool   bigint_sign(ConstBigInt z);
const BigInt_DataBlock* bigint_data(ConstBigInt z);

BigInt bigint_abs(ConstBigInt z, BigInt* out);
BigInt bigint_neg(ConstBigInt z, BigInt* out);

BigInt bigint_create_zero();
BigInt bigint_create_small(SmallInt v);
BigInt bigint_create_usmall(USmallInt v);

BigInt bigint_set_zero(BigInt* z);
BigInt bigint_set_small(BigInt* z, SmallInt v);
BigInt bigint_set_usmall(BigInt* z, USmallInt v);

BigInt bigint_copy(BigInt* dst, ConstBigInt src);

int bigint_ucmp(ConstBigInt a, ConstBigInt b);
int bigint_cmp(ConstBigInt a, ConstBigInt b);

// returns the first "big digit" without a sign, as a small unsigned integer
USmallInt bigint_usmall(ConstBigInt a);

// returns the first "big digit" with a sign, as a small signed integer
SmallInt bigint_small(ConstBigInt a);

// unsigned comparison of big integer with small integer
int bigint_ucmp_small(ConstBigInt a, USmallInt b);

// signed comparison of big integer with small integer
int bigint_cmp_small(ConstBigInt a, SmallInt b);

// unsigned addition
BigInt bigint_uadd(ConstBigInt a, ConstBigInt b, BigInt* out);

// unsigned subtraction
BigInt bigint_usub(ConstBigInt a, ConstBigInt b, BigInt* out);

// unsigned multiplication
BigInt bigint_umul(ConstBigInt a, ConstBigInt b, BigInt* out);

// unsigned division
BigInt bigint_udiv(ConstBigInt a, ConstBigInt b, BigInt* out, BigInt* rem);

// signed addition
BigInt bigint_add(ConstBigInt a, ConstBigInt b, BigInt* out);

// signed subtraction
BigInt bigint_sub(ConstBigInt a, ConstBigInt b, BigInt* out);

// signed multiplication
BigInt bigint_mul(ConstBigInt a, ConstBigInt b, BigInt* out);

// signed division
BigInt bigint_div(ConstBigInt a, ConstBigInt b, BigInt* out, BigInt* rem);

BigInt bigint_lshift(ConstBigInt z, size_t shift, BigInt* out);
BigInt bigint_rshift(ConstBigInt z, size_t shift, BigInt* out);
BigInt bigint_lshift_blocks(ConstBigInt z, size_t shift, BigInt* out);
BigInt bigint_rshift_blocks(ConstBigInt z, size_t shift, BigInt* out);

char* bigint_str(ConstBigInt z, BigInt_FormatSpec bifs);
char* bigint_hex_str(ConstBigInt z, BigInt_FormatSpec bifs);
char* bigint_dec_str(ConstBigInt z, BigInt_FormatSpec bifs);

BigInt bigint_scan(const char* str, BigInt* out);

int bigint_fprintf(FILE* file, const char* const format, ...);
int bigint_printf(const char* const format, ...);
