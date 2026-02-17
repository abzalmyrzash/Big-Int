#pragma once
#include "bigint_alias.h"

static inline size_t size_without_zeros(const Block* a, size_t size) {
	while (size > 0) {
		size--;
		if (a[size] != 0) return size + 1;
	}
	return size;
}

