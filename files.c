#include "files.h"
#include <stdio.h>
#include <windows.h>

bool printFile(const char* filename) {
	FILE* file = fopen(filename, "r");
	if (!file) {
		return false;
	}
	char c;
	while((c = fgetc(file)) != EOF) {
		putchar(c);
	}
	fclose(file);
	return true;
}

char* getBasename(char* path) {
	char* ptr = path + strlen(path);
	while (ptr != path) {
		if (*ptr == '/' || *ptr == '\\') {
			return ptr + 1;
		}
		ptr--;
	}
	return ptr;
}

char* getExtension(char* filename) {
	char* ptr = filename + strlen(filename);
	while (ptr != filename) {
		if (*ptr == '.') {
			return ptr;
		}
		ptr--;
	}
	return NULL;
}

int mkdir(const char* dirname) {
	if (!CreateDirectory(dirname, NULL)) {
		if(ERROR_ALREADY_EXISTS == GetLastError())
			return DIR_ALREADY_EXISTS;
		return DIR_ERROR;
	}
	return DIR_SUCCESS;
}

