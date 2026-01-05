#include "bigint.h"
#include "bigint_rand.h"
#include "fib.h"
#include "factorial.h"
#include <assert.h>
#include <stdio.h>
#include <stdbool.h>
#include <time.h>

int run(BigInt (*f)(unsigned int, BigInt*));

int main()
{
	/*
	run(bigFib);
	return 0;
	bigint_init();
	BigInt A = NULL;
	BigInt B = NULL;
	bigint_rand(&A, 10'000);
	bigint_rand(&B, 10'000);
	printf("%zu %zu\n", bigint_size(A), bigint_size(B));
	BigInt C = bigint_alloc(bigint_size(A) + bigint_size(B));
	bigint_memstat();
	clock_t start = clock();
	for (int i = 0; i < 10000; i++) {
		bigint_mul(A, B, &C);
	}
	float time = (float)(clock() - start) / CLOCKS_PER_SEC;
	printf("Time: %f seconds\n", time);
	bigint_free(A);
	bigint_free(B);
	bigint_free(C);
	bigint_memstat();
	return 0;
	*/

	bigint_init();
	bigint_finish();
	return 0;
	size_t p;
	scanf("%zu", &p);
	set_newton_precision(p);
	return 0;

	BigInt A = NULL;
	BigInt B = NULL;
	BigInt C = NULL;
	BigInt res = NULL;
	BigInt rem = NULL;

	constexpr int buffer_size = 102400;
	char buffer[buffer_size];

	printf("Enter a: ");
	fgets(buffer, buffer_size, stdin);
	bigint_scan(buffer, &A);
	bigint_printf("%p0_x\n", A);
	// bigint_printf("%0_d\n", A);
	printf("Size: %zu\n", bigint_size(A));

	printf("Enter b: ");
	fgets(buffer, buffer_size, stdin);
	bigint_scan(buffer, &B);
	bigint_printf("%p0_x\n", B, B);
	// bigint_printf("%0_d\n", B, B);
	printf("Size: %zu\n", bigint_size(B));

	printf("Enter c: ");
	fgets(buffer, buffer_size, stdin);
	bigint_scan(buffer, &C);
	bigint_printf("%p0_x\n", C, C);
	// bigint_printf("%0_d\n", B, B);
	printf("Size: %zu\n", bigint_size(C));

	bigint_powmod(A, B, C, &res);
	printf("Size: %zu\n", bigint_size(res));
	fflush(stdout);
	bigint_printf("%p0_x\n", res);
	fflush(stdout);
	bigint_printf("%d\n", res);

	return 0;
	res = bigint_alloc(bigint_size(A) + bigint_size(B));

	bigint_add(A, B, &res);
	printf("A + B = ");
	bigint_printf("%d\n", res);

	assert(bigint_cmp(A, bigint_sub(res, B, &res)) == 0 || ({
				bigint_printf("%d\n", res);
				}));

	bigint_sub(A, B, &res);
	printf("A - B = ");
	bigint_printf("%d\n", res);

	assert(bigint_cmp(A, bigint_add(res, B, &res)) == 0);

	bigint_mul(A, B, &res);
	printf("A * B = ");
	bigint_printf("%d\n", res);

	assert(bigint_cmp(A, bigint_div(res, B, &res, &rem)) == 0);
	assert(bigint_cmp_small(rem, 0) == 0);

	bigint_div(A, B, &res, &rem);
	printf("A / B = ");
	bigint_printf("%d\n", res);
	printf("A %% B = ");
	bigint_printf("%d\n", rem);

	assert(bigint_ucmp(B, rem) > 0);

	bigint_add(rem, bigint_mul(res, B, &res), &res);
	bigint_printf("%d\n", res);

	assert(bigint_cmp(res, A) == 0);

	bigint_free(A);
	bigint_free(B);
	bigint_free(res);
	bigint_free(rem);
	bigint_memstat();
	bigint_finish();
	return 0;
}

int run(BigInt (*f)(unsigned int, BigInt*)) {
	bigint_init();
	BigInt res = NULL;
	while(1) {
		int n;
		scanf("%d", &n);
		if (n < 0) break;

		clock_t start = clock();
		f(n, &res);
		double time = (float)(clock() - start) / CLOCKS_PER_SEC;

		start = clock();
		bigint_printf("%px\n\n", res);
		double time2 = (float)(clock() - start) / CLOCKS_PER_SEC;

		printf("Size: %zu\n", bigint_size(res));
		printf("Cap: %zu\n", bigint_cap(res));
		printf("Time to calculate: %f seconds\n", time);
		printf("Time to print hexadecimal: %f seconds\n", time2);
		//bigint_memstat();

		bool convert_to_decimal = true;
		if (bigint_size(res) > 200) {
			printf("Convert to decimal? (y/n)");
			char c;
			scanf(" %c", &c);
			convert_to_decimal = (c == 'y' || c == 'Y');
		}

		if (convert_to_decimal) {
			start = clock();
			bigint_printf("\n%d\n\n", res);
			double time3 = (float)(clock() - start) / CLOCKS_PER_SEC;
			printf("Time to convert to decimal: %f seconds\n", time3);
			//bigint_memstat();
		}
		printf("\n");
	}
	bigint_free(res);
	bigint_memstat();
	bigint_finish();
	return 0;
}
