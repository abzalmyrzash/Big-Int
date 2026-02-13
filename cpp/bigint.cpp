#include "bigint.hpp"

extern "C" {
	#include "..\bigint.h"
};

static int cnt = 0;
bool BigIntClass::is_log = false;

static BigInt ONE;

void BigIntClass::init() {
	bigint_init();
	ONE = bigint_create_usmall(1);
}

void BigIntClass::finish() {
	bigint_finish();
	bigint_free(ONE);
}

void BigIntClass::log(bool y) {
	bigint_log_memory(y);
	is_log = y;
}

void BigIntClass::reserve(size_t cap, bool keep_data) {
	bigint_reserve(&x, cap, keep_data);
}

void BigIntClass::destroy() {
	bigint_free(x);
	x = NULL;
}

BigIntClass::BigIntClass() {
	id = cnt++;
	if (is_log) std::cout << "Default constructor " << id << std::endl;
	x = NULL;
}

BigIntClass::BigIntClass(SmallInt v) {
	id = cnt++;
	if (is_log) std::cout << "SmallInt constructor " << id << std::endl;
	x = bigint_create_small(v);
}

BigIntClass::BigIntClass(const BigIntClass& b) {
	id = cnt++;
	if (is_log) std::cout << "copy constructor " << id << " = " << b.id << std::endl;
	x = NULL;
	bigint_copy(&x, b.x);
}

BigIntClass::BigIntClass(BigIntClass&& b) {
	id = cnt++;
	if (is_log) std::cout << "move constructor " << id << " = " << b.id << std::endl;
	x = b.x;
	b.x = NULL;
}

BigIntClass::~BigIntClass() {
	if (is_log) std::cout << "destructor " << id << std::endl;
	bigint_free(x);
}

BigIntClass& BigIntClass::operator=(const BigIntClass& b) {
	if (is_log) std::cout << "copy assignment " << id << " = " << b.id << std::endl;
	bigint_copy(&x, b.x);
	return *this;
}

BigIntClass& BigIntClass::operator=(BigIntClass&& b) {
	if (is_log) std::cout << "move assignment " << id << " = " << b.id << std::endl;
	if (x != b.x) {
		bigint_free(x);
	}
	x = b.x;
	b.x = NULL;
	return *this;
}

BigIntClass& BigIntClass::operator=(SmallInt v) {
	bigint_set_small(&x, v);
	return *this;
}

BigIntClass& BigIntClass::operator=(USmallInt v) {
	bigint_set_usmall(&x, v);
	return *this;
}

BigIntClass& BigIntClass::operator=(signed v) {
	bigint_set_small(&x, v);
	return *this;
}

BigIntClass& BigIntClass::operator=(unsigned v) {
	bigint_set_usmall(&x, v);
	return *this;
}

bool BigIntClass::operator==(const BigIntClass& b) {
	return bigint_cmp(x, b.x) == 0;
}

bool BigIntClass::operator>(const BigIntClass& b) {
	return bigint_cmp(x, b.x) > 0;
}

bool BigIntClass::operator<(const BigIntClass& b) {
	return bigint_cmp(x, b.x) < 0;
}

bool BigIntClass::operator>=(const BigIntClass& b) {
	return bigint_cmp(x, b.x) >= 0;
}

bool BigIntClass::operator<=(const BigIntClass& b) {
	return bigint_cmp(x, b.x) <= 0;
}

bool BigIntClass::operator==(SmallInt v) {
	return bigint_cmp_small(x, v) == 0;
}

bool BigIntClass::operator>(SmallInt v) {
	return bigint_cmp_small(x, v) > 0;
}

bool BigIntClass::operator<(SmallInt v) {
	return bigint_cmp_small(x, v) < 0;
}

bool BigIntClass::operator>=(SmallInt v) {
	return bigint_cmp_small(x, v) >= 0;
}

bool BigIntClass::operator<=(SmallInt v) {
	return bigint_cmp_small(x, v) <= 0;
}

BigIntClass BigIntClass::operator+(const BigIntClass& b) {
	if (is_log) std::cout << "add" << std::endl;
	BigIntClass c;
	bigint_add(x, b.x, &c.x);
	return c;
}

BigIntClass BigIntClass::operator-(const BigIntClass& b) {
	if (is_log) std::cout << "sub" << std::endl;
	BigIntClass c;
	bigint_sub(x, b.x, &c.x);
	return c;
}

BigIntClass BigIntClass::operator*(const BigIntClass& b) {
	if (is_log) std::cout << "mul" << std::endl;
	BigIntClass c;
	bigint_mul(x, b.x, &c.x);
	return c;
}

BigIntClass BigIntClass::operator/(const BigIntClass& b) {
	if (is_log) std::cout << "div" << std::endl;
	BigIntClass c;
	bigint_div(x, b.x, &c.x, NULL);
	return c;
}

BigIntClass BigIntClass::operator%(const BigIntClass& b) {
	if (is_log) std::cout << "mod" << std::endl;
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
	if (BigIntClass::is_log) std::cout << "read" << std::endl;
	std::string s;
	ins >> s;
	bigint_scan(s.data(), s.size(), &z.x);
	return ins;
}

BigIntClass& BigIntClass::operator+=(const BigIntClass& b) {
	bigint_add(x, b.x, &x);
	return *this;
}

BigIntClass& BigIntClass::operator-=(const BigIntClass& b) {
	bigint_sub(x, b.x, &x);
	return *this;
}

BigIntClass& BigIntClass::operator*=(const BigIntClass& b) {
	bigint_mul(x, b.x, &x);
	return *this;
}

BigIntClass& BigIntClass::operator/=(const BigIntClass& b) {
	bigint_div(x, b.x, &x, NULL);
	return *this;
}

BigIntClass& BigIntClass::operator%=(const BigIntClass& b) {
	bigint_div(x, b.x, NULL, &x);
	return *this;
}

BigIntClass& BigIntClass::add(const BigIntClass& a, const BigIntClass& b) {
	bigint_add(a.x, b.x, &x);
	return *this;
}

BigIntClass& BigIntClass::sub(const BigIntClass& a, const BigIntClass& b) {
	bigint_sub(a.x, b.x, &x);
	return *this;
}

BigIntClass& BigIntClass::mul(const BigIntClass& a, const BigIntClass& b) {
	bigint_mul(a.x, b.x, &x);
	return *this;
}

BigIntClass& BigIntClass::div(const BigIntClass& a, const BigIntClass& b) {
	bigint_div(a.x, b.x, &x, NULL);
	return *this;
}

BigIntClass& BigIntClass::mod(const BigIntClass& a, const BigIntClass& b) {
	bigint_div(a.x, b.x, NULL, &x);
	return *this;
}

BigIntClass& BigIntClass::operator++() {
	bigint_add(x, ONE, &x);
	return *this;
}

BigIntClass  BigIntClass::operator++(int) {
	BigIntClass copy(*this);
	++(*this);
	return copy;
}

BigIntClass& BigIntClass::operator--() {
	bigint_sub(x, ONE, &x);
	return *this;
}

BigIntClass  BigIntClass::operator--(int) {
	BigIntClass copy(*this);
	--(*this);
	return copy;
}

