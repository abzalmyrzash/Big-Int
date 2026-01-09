CC_FLAGS := -g -O3 --std=c2x -lm
CC := gcc

build/main: main.c build/fib.o build/factorial.o build/bigint.o build/bigint_rand.o build/arena.o utils.h
	$(CC) $(CC_FLAGS) main.c build/*.o -o build/main

build/bigint.o: bigint.c bigint.h bigint_impl.h bigint_impl_basic.h bigint_params.h utils.h
	$(CC) $(CC_FLAGS) -c bigint.c -o build/bigint.o

build/fib.o: fib.c fib.h utils.h
	$(CC) $(CC_FLAGS) -c fib.c -o build/fib.o

build/factorial.o: factorial.c factorial.h utils.h
	$(CC) $(CC_FLAGS) -c factorial.c -o build/factorial.o

build/bigint_rand.o: bigint_rand.c bigint_rand.h utils.h
	$(CC) $(CC_FLAGS) -c bigint_rand.c -o build/bigint_rand.o

build/arena.o: arena.c arena.h utils.h
	$(CC) $(CC_FLAGS) -c arena.c -o build/arena.o

clean:
	del build\*
