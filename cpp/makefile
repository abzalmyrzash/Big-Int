debug:
	clang++ -std=c++23 -g main.cpp bigint.cpp rsa.cpp base64.cpp ..\debug\bigint.o ..\debug\arena.o ..\debug\bigint_rand.o ..\debug\bigint_crypt.o ..\debug\clipboard.o ..\debug\prompt.o ..\debug\files.o -o main.exe -lws2_32 -luser32

release:
	clang++ -DNDEBUG -std=c++23 -g main.cpp bigint.cpp rsa.cpp base64.cpp ..\release\bigint.o ..\release\arena.o ..\release\bigint_rand.o ..\release\bigint_crypt.o ..\release\clipboard.o ..\release\prompt.o ..\release\files.o -o main.exe -lws2_32 -luser32
