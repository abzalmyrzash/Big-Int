CC_FLAGS := -O0 -g --std=c23
CC := gcc

build/main: main.c build/fib.o build/factorial.o build/bigint.o
	$(CC) $(CC_FLAGS) main.c build/\*.o -o build/main

build/fib.o: fib.c fib.h build\bigint.o
	$(CC) $(CC_FLAGS) -c fib.c build/bigint.o -o build/fib.o

build/factorial.o: factorial.c factorial.h build\bigint.o
	$(CC) $(CC_FLAGS) -c factorial.c build/bigint.o -o build/factorial.o

build/bigint.o: bigint.c bigint.h bigint_params.h
	$(CC) $(CC_FLAGS) -c bigint.c -o build/bigint.o
