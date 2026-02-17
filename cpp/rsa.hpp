#pragma once
#include "bigint.hpp"
#include <vector>

void generateKeys(const size_t width, BigIntClass& n, BigIntClass& d, BigIntClass& e, bool showVal = false);

void divide(const std::string& str, size_t bitsPerBlock,
		std::vector<BigIntClass>& out);

void divide(const char* str, const size_t strLen, size_t bitsPerBlock,
		std::vector<BigIntClass>& out);

void crypt(const std::vector<BigIntClass>& blocks, BigIntClass n, BigIntClass e,
		std::vector<BigIntClass>& out);

void merge(const std::vector<BigIntClass>& blocks, size_t bitsPerBlock,
		std::string& out);

bool testKeys(BigIntClass n, BigIntClass d, BigIntClass e);

