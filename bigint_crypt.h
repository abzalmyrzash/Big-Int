#pragma once
#include "bigint.h"

bool crypt_init();
BigInt bigint_crypt_rand(size_t width, BigInt* out_ptr);
bool crypt_finish();

