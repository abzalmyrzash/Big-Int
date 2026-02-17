#pragma comment(lib, "advapi32.lib")
#include <windows.h>
#include "bigint_crypt.h"
#include "bigint_impl_basic.h"

static HCRYPTPROV hCryptProv;

bool crypt_init() {
	return CryptAcquireContext(
			&hCryptProv,
			NULL,
			"Microsoft Base Cryptographic Provider v1.0",
			PROV_RSA_FULL,
			CRYPT_VERIFYCONTEXT);
}

BigInt bigint_crypt_rand(size_t width, BigInt* out_ptr) {
	size_t blocks = width / BLOCK_WIDTH;
	size_t rem_bits = width % BLOCK_WIDTH;
	size_t cap = blocks + (rem_bits > 0);
	BigInt out = bigint_reserve(out_ptr, cap, KEEP_DATA);
	if (!out) return NULL;

	size_t bytes = cap * sizeof(Block);
	if(CryptGenRandom(
	   hCryptProv,
	   bytes, 
	   (BYTE*)out->data))
	{
		if (rem_bits) out->data[cap - 1] &= ((1ULL << rem_bits) - 1);
		out->size = size_without_zeros(out->data, cap);
		return out;
	}
	else
	{
		 printf("Error during CryptGenRandom.\n");
		 return NULL;
	}
}

bool crypt_finish() {
	return CryptReleaseContext(hCryptProv, 0);
}

