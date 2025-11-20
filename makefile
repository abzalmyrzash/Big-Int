build\main: *.c *.h makefile
	gcc -O3 --std=c23 *.c -o build\main
