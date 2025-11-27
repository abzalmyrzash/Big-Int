CC := cl

build/main: main.c build\fib.o build\factorial.o build\bigint.o
	$(CC) $(CC_FLAGS) main.c build/\*.o /link /out:build\main

build/fib.o: fib.c fib.h build\bigint.o
	$(CC) $(CC_FLAGS) /c fib.c build/bigint.o /link /out:build\fib.o

build/factorial.o: factorial.c factorial.h build\bigint.o
	$(CC) $(CC_FLAGS) /c factorial.c build/bigint.o /link /out:build\factorial.o

build/bigint.o: bigint.c bigint.h bigint_params.h
	$(CC) $(CC_FLAGS) /c bigint.c /out:build/bigint.o
