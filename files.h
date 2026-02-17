#pragma once
#include <stdbool.h>

bool printFile(const char* filename);

char* getBasename(char* path);

char* getExtension(char* filename);

#define DIR_ALREADY_EXISTS  1
#define DIR_SUCCESS         0
#define DIR_ERROR          -1

int mkdir(const char* dirName);
