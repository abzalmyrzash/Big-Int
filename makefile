CC_FLAGS := -g -Wall --std=c23
CC := gcc

build/main: main.c build/fib.o build/factorial.o build/bigint.o build/bigint_rand.o makefile
	$(CC) $(CC_FLAGS) main.c build/\*.o -o build/main

build/fib.o: fib.c fib.h build/bigint.o makefile
	$(CC) $(CC_FLAGS) -c fib.c build/bigint.o -o build/fib.o

build/factorial.o: factorial.c factorial.h build/bigint.o makefile
	$(CC) $(CC_FLAGS) -c factorial.c build/bigint.o -o build/factorial.o

build/bigint_rand.o: bigint_rand.c bigint_rand.h build/bigint.o makefile
	$(CC) $(CC_FLAGS) -c bigint_rand.c build/bigint.o -o build/bigint_rand.o

build/bigint.o: bigint.c bigint.h bigint_impl.h bigint_params.h makefile
	$(CC) $(CC_FLAGS) -c bigint.c -o build/bigint.o

clean:
	del build\*
