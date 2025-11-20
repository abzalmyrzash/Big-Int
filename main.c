#include "bigint.h"
#include "fib.h"
#include "factorial.h"
#include <assert.h>
#include <stdio.h>
#include <stdbool.h>
#include <time.h>

int main()
{
	/*
	BigInt A = NULL;
	BigInt B = NULL;
	BigInt res = NULL;
	BigInt rem = NULL;

	printf("Enter a: ");
	constexpr int buffer_size = 256;
	char buffer[buffer_size];
	fgets(buffer, buffer_size, stdin);
	bigint_scan(buffer, &A);

	bigint_printf("%p0_x\n%0_d\n", A, A);

	printf("Enter b: ");
	fgets(buffer, buffer_size, stdin);
	bigint_scan(buffer, &B);
	bigint_printf("%p0_x\n%0_d\n", B, B);

	bigint_add(A, B, &res);
	printf("A + B = ");
	bigint_printf("%d\n", res);

	assert(bigint_cmp(A, bigint_sub(res, B, &res)) == 0);

	bigint_sub(A, B, &res);
	printf("A - B = ");
	bigint_printf("%d\n", res);

	assert(bigint_cmp(A, bigint_add(res, B, &res)) == 0);

	bigint_mul(A, B, &res);
	printf("A * B = ");
	bigint_printf("%d\n", res);

	assert(bigint_cmp(A, bigint_div(res, B, &res, NULL)) == 0);

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
	*/

	BigInt res = bigint_alloc(10000);
	while(1) {
		int n;
		scanf("%d", &n);
		if (n < 0) return 0;

		clock_t start = clock();
		bigFib(n, &res);
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

	return 0;
}
