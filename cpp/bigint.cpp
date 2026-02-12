#include "bigint.hpp"

extern "C" {
	#include "..\bigint.h"
};

BigIntClass::BigIntClass() {
	x = NULL;
}

BigIntClass::BigIntClass(SmallInt v) {
	x = bigint_create_small(v);
}

BigIntClass::~BigIntClass() {
	bigint_free(x);
}

BigIntClass BigIntClass::operator=(SmallInt v) {
	bigint_set_small(&x, v);
	return *this;
}

BigIntClass BigIntClass::operator+(const BigIntClass& b) {
	BigIntClass c;
	bigint_add(x, b.x, &c.x);
	return c;
}

BigIntClass BigIntClass::operator-(const BigIntClass& b) {
	BigIntClass c;
	bigint_sub(x, b.x, &c.x);
	return c;
}

BigIntClass BigIntClass::operator*(const BigIntClass& b) {
	BigIntClass c;
	bigint_mul(x, b.x, &c.x);
	return c;
}

BigIntClass BigIntClass::operator/(const BigIntClass& b) {
	BigIntClass c;
	bigint_div(x, b.x, &c.x, NULL);
	return c;
}

BigIntClass BigIntClass::operator%(const BigIntClass& b) {
	BigIntClass c = BigIntClass(0);
	bigint_div(x, b.x, NULL, &c.x);
	return c;
}

std::ostream& operator << (std::ostream& outs, const BigIntClass& z) {
	BigInt_FormatSpec bifs = { 0 };
	bifs.base = 10;
	char* str = bigint_dec_str(z.x, bifs);
	outs << str;
	free(str);
	return outs;
}

std::istream& operator >> (std::istream& ins, BigIntClass& z) {
	std::string s;
	ins >> s;
	bigint_scan(s.data(), s.size(), &z.x);
	return ins;
}

