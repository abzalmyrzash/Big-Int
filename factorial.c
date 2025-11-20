#include "factorial.h"

BigInt bigFactorial(unsigned int n, BigInt* out_ptr) {
	BigInt out = bigint_set_usmall(out_ptr, 1);
	if (n <= 1) return out;
	for (int i = 2; i <= n; i++) {
		BigInt j = bigint_create_usmall(i);
		out = bigint_umul(out, j, out_ptr);
		bigint_free(j);
	}
	return out;
}
