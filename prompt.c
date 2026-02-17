#include "prompt.h"
#include <stdio.h>
#include <ctype.h>
#include <conio.h>

char getOption() {
	char c = toupper(getchar());
	if (c != '\n') while (getchar() != '\n');
	return c;
}

void waitChar(char c, const char* msg) {
	printf("%s", msg);
	if (c == ANY_KEY) _getch();
	else while(_getch() != c);
}

void pressAnyKeyToContinue() {
	waitChar(ANY_KEY, "Press any key to continue...");
	printf("\n");
}

