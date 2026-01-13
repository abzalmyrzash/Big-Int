#include "bigint.h"
#include "bigint_rand.h"
#include "fib.h"
#include "factorial.h"
#include "utils.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>

int run(BigInt (*f)(unsigned int, BigInt*));

int main()
{
	clock_t start;
	float elapsed_time;
	size_t size;

	bigint_init();

	BigInt A = NULL;
	BigInt B = NULL;
	BigInt C = NULL;
	BigInt res = NULL;
	BigInt rem = NULL;
	BigInt recpr = NULL;

	/*
	for (int i = 0; i < 1000000; i++) {
		size_t a_width = rand();
		size_t b_width = rand();
		bigint_rand(&A, a_width);
		bigint_rand(&B, b_width);
		// printf("%zu %zu\n", bigint_size(A), bigint_size(B));
		start = clock();
		bigint_mul(A, B, &C);
		elapsed_time = (float)(clock() - start) / CLOCKS_PER_SEC;
		// printf("Time: %f seconds\n", elapsed_time);
	}
	bigint_free(A);
	bigint_free(B);
	bigint_free(C);
	bigint_memstat();
	return 0;
	*/

	BigInt num = NULL;

	/*
	FILE* file = fopen("dec.txt", "r");
	char str[100000];
	fscanf(file, "%s", str);
	start = clock();
	bigint_scan(str, &num);
	elapsed_time = (float)(clock() - start) / CLOCKS_PER_SEC;
	fclose(file);
	printf("Time: %.3f s\n", elapsed_time);

	start = clock();
	bigint_printf("%d\n", num);
	elapsed_time = (float)(clock() - start) / CLOCKS_PER_SEC;
	printf("Time: %.3f s\n", elapsed_time);
	*/

	size_t max_width = 100'000;
	// scanf("%zu", &max_width);
	bigint_rand(&num, max_width);

	FILE* file = fopen("dec_spaces.txt", "w");
	start = clock();
	bigint_fprintf(file, "%_d\n", num);
	elapsed_time = (float)(clock() - start) / CLOCKS_PER_SEC;
	fclose(file);
	printf("Time: %.3f s\n", elapsed_time);
	bigint_finish();

	return 0;

	scanf("%zu", &size);
	size_t precision = size * BIGINT_BLOCK_WIDTH;
	bigint_rand(&A, precision);
	bigint_rand(&B, precision);

	start = clock();
	bigint_mul(A, B, &res);
	elapsed_time = (float)(clock() - start) / CLOCKS_PER_SEC;
	printf("Time: %.3f s\n", elapsed_time);
	printf("Size: %zu\n", bigint_size(res));
	return 0;

	/*
	constexpr int buffer_size = 102400;
	char buffer[buffer_size];

	printf("Enter a: ");
	fgets(buffer, buffer_size, stdin);
	bigint_scan(buffer, &A);

	printf("Enter b: ");
	fgets(buffer, buffer_size, stdin);
	bigint_scan(buffer, &B);
	*/


	bigint_printf("%p0_x\n", A);
	bigint_printf("%p0_x\n", B);
	bigint_printf("%0_d\n", A);
	printf("Size: %zu\n", bigint_size(A));
	bigint_printf("%0_d\n", B, B);
	printf("Size: %zu\n", bigint_size(B));

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

	bigint_recpr(B, precision, &recpr);
	bigint_printf("1 / B = ");
	bigint_printf("%px\n", recpr);
	bigint_div_recpr(A, B, recpr, precision, &res, &rem);
	bigint_printf("A / B (using reciprocal) = %d\n", res);
	bigint_printf("A %% B (using reciprocal) = %d\n", rem);

	bigint_free(A);
	bigint_free(B);
	bigint_free(res);
	bigint_free(rem);
	bigint_free(recpr);
	bigint_memstat();
	bigint_finish();
	return 0;
}

int run(BigInt (*f)(unsigned int, BigInt*)) {
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
	return 0;
}
