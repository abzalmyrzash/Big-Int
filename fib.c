#include "fib.h"
#include <assert.h>
#include <math.h>
#include "elog.h"

uint64_t fib(unsigned int n)
{
	if (n == 0) return 0;
	if (n == 1) return 1;
	uint64_t a = 0;
	uint64_t b = 1;
	uint64_t tmp;
	for (unsigned int i = 2; i <= n; i++) {
		tmp = b;
		b += a;
		a = tmp;
	}
	return b;
}

//constexpr double GOLDEN_RATIO = 1.61803398875;
constexpr double LOG2_GOLDEN_RATIO = 0.69424191363; 

BigInt bigFib(unsigned int n, BigInt* out_ptr)
{
	if (n == 0) {
		return bigint_set_zero(out_ptr);
	}
	if (n == 1 || n == 2) {
		return bigint_set_usmall(out_ptr, 1);
	}
	BigInt a, out, tmp;

	const size_t cap = ceil(n * LOG2_GOLDEN_RATIO / BIGINT_BLOCK_WIDTH);

	a = bigint_alloc(cap);
	a = bigint_set_usmall(&a, 1);
	out = bigint_reserve(out_ptr, cap, false);
	out = bigint_set_usmall(out_ptr, 1);
	tmp = bigint_alloc(cap);

	for (unsigned int i = 3; i <= n; i++) {
		bigint_copy(&tmp, out);
		out = bigint_uadd(out, a, out_ptr);
		if (!out) {
			ELOG_STR("failed to add");
			return NULL;
		}
		bigint_copy(&a, tmp);
	}

	bigint_free(a);
	bigint_free(tmp);
	return out;
}
