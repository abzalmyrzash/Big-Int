#pragma once
#include "bigint.hpp"

inline BigIntClass gcd(BigIntClass a, BigIntClass b) {
	while (b != 0) {
		BigIntClass tmp = b;
		b.mod(a, b);
		a = std::move(tmp);
	}
	return a;
}

inline BigIntClass lcm(const BigIntClass& a, const BigIntClass& b) {
	return a * b / gcd(a, b);
}

/*
inline BigIntClass powmod(const BigIntClass& M, BigIntClass e, BigIntClass n) {
	BigIntClass C = 1;
	const size_t bits = e.width();
	for (size_t i = bits - 1; ; i--) {
		C.mod((C * C), n);
		if (e.bit(i)) {
			C.mod((C * M), n);
		}
		if (i == 0) break;
	}
	return C;
}
*/

inline BigIntClass powmod(const BigIntClass& _M, const BigIntClass& e, const BigIntClass& n) {
	BigIntClass M = _M % n;
	BigIntClass C = 1;
	const size_t bits = e.width();
	const size_t prec = 2 * n.width();
	BigIntClass recpr = n.recpr(prec);
	for (size_t i = bits - 1; ; i--) {
		C.mod((C * C), n, recpr, prec);
		if (e.bit(i)) {
			C.mod((C * M), n, recpr, prec);
		}
		if (i == 0) break;
	}
	return C;
}

