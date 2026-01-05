#pragma once
#include "bigint_alias.h"

static inline size_t size_without_zeros(const Block* a, size_t size) {
	while (size > 0) {
		if (a[size - 1] != 0) break;
		size--;
	}
	return size;
}

