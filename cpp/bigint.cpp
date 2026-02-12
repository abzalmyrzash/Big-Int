#include "bigint.hpp"

extern "C" {
	#include "..\bigint.h"
};

static int cnt = 0;

BigIntClass::BigIntClass() {
	x = NULL;
	id = cnt++;
}

BigIntClass::BigIntClass(SmallInt v) {
	x = bigint_create_small(v);
	id = cnt++;
}

BigIntClass::BigIntClass(const BigIntClass& b) {
	x = NULL;
	bigint_copy(&x, b.x);
	id = cnt++;
}

BigIntClass::BigIntClass(BigIntClass&& b) {
	x = b.x;
	b.x = NULL;
	id = cnt++;
}

BigIntClass::~BigIntClass() {
	// std::cout << "destructor " << id << std::endl;
	bigint_free(x);
}

BigIntClass BigIntClass::operator=(const BigIntClass& b) {
	bigint_copy(&x, b.x);
	// std::cout << "copy" << std::endl;
	return *this;
}

BigIntClass BigIntClass::operator=(BigIntClass&& b) {
	if (x != b.x) {
		bigint_free(x);
	}
	x = b.x;
	b.x = NULL;
	// std::cout << "move " << id << " = " << b.id << std::endl;
	return *this;
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
	BigIntClass c;
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

