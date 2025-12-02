#include "factorial.h"

BigInt bigFactorial(unsigned int n, BigInt* out_ptr) {
	BigInt out = bigint_set_usmall(out_ptr, 1);
	if (n <= 1) return out;
	BigInt j = NULL;
	BigInt tmp = NULL;
	for (int i = 2; i <= n; i++) {
		bigint_set_usmall(&j, i);
		bigint_copy(&tmp, out);
		out = bigint_umul(tmp, j, out_ptr);
	}
	bigint_free(j);
	bigint_free(tmp);
	return out;
}
