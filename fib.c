#include "fib.h"
#include <assert.h>

uint64_t fib(unsigned int n)
{
	if (n == 0) return 0;
	if (n == 1) return 1;
	uint64_t a = 0;
	uint64_t b = 1;
	uint64_t tmp;
	for (int i = 2; i <= n; i++) {
		tmp = b;
		b += a;
		a = tmp;
	}
	return b;
}

BigInt bigFib(unsigned int n, BigInt* out_ptr)
{
	if (n == 0) {
		return bigint_set_zero(out_ptr);
	}
	if (n == 1 || n == 2) {
		return bigint_set_usmall(out_ptr, 1);
	}
	BigInt a = bigint_create_usmall(1);
	BigInt out = bigint_set_usmall(out_ptr, 1);
	BigInt tmp = NULL;
	for (int i = 3; i <= n; i++) {
		bigint_copy(&tmp, out);
		out = bigint_uadd(out, a, out_ptr);
		bigint_copy(&a, tmp);
	}
	bigint_free(a);
	bigint_free(tmp);
	return out;
}
