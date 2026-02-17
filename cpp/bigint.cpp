#include "bigint.hpp"
#include "..\utils.h"

extern "C" {
	#include "..\bigint.h"
	#include "..\bigint_rand.h"
	#include "..\bigint_crypt.h"
};

#include <cassert>

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

size_t BigIntClass::size() const {
	return bigint_size(x);
}

size_t BigIntClass::width() const {
	return bigint_width(x);
}

size_t BigIntClass::cap() const {
	return bigint_cap(x);
}

size_t width_to_size(size_t width) {
	return ceil_div(width, BigIntClass::BLOCK_WIDTH);
}

using Block = BigIntClass::Block;

const Block* BigIntClass::data() const {
	return bigint_data(x);
}

void BigIntClass::reserve(size_t cap, bool keep_data) {
	bigint_reserve(&x, cap, keep_data);
	assert(bigint_is_valid(x));
}

void BigIntClass::destroy() {
	bigint_free(x);
	x = NULL;
}

BigIntClass::BigIntClass() {
	x = NULL;
}

BigIntClass::BigIntClass(SmallInt v) {
	x = bigint_create_small(v);
	assert(bigint_is_valid(x));
}

BigIntClass::BigIntClass(USmallInt v) {
	x = bigint_create_usmall(v);
	assert(bigint_is_valid(x));
}

BigIntClass::BigIntClass(signed v) {
	x = bigint_create_small(v);
	assert(bigint_is_valid(x));
}

BigIntClass::BigIntClass(unsigned v) {
	x = bigint_create_usmall(v);
	assert(bigint_is_valid(x));
}

BigIntClass::BigIntClass(const BigIntClass& b) {
	x = NULL;
	bigint_copy(&x, b.x);
	assert(bigint_is_valid(x));
	if (is_log) std::cout << "copy constructor" << std::endl;
}

BigIntClass::BigIntClass(BigIntClass&& b) {
	x = b.x;
	b.x = NULL;
	assert(bigint_is_valid(x));
	if (is_log) std::cout << "move constructor" << std::endl;
}

BigIntClass::~BigIntClass() {
	bigint_free(x);
}

BigIntClass& BigIntClass::operator=(const BigIntClass& b) {
	bigint_copy(&x, b.x);
	assert(bigint_is_valid(x));
	if (is_log) std::cout << "copy assignment" << std::endl;
	return *this;
}

BigIntClass& BigIntClass::operator=(BigIntClass&& b) {
	if (x != b.x) {
		bigint_free(x);
	}
	x = b.x;
	b.x = NULL;
	assert(bigint_is_valid(x));
	if (is_log) std::cout << "move assignment" << std::endl;
	return *this;
}

BigIntClass& BigIntClass::operator=(SmallInt v) {
	bigint_set_small(&x, v);
	assert(bigint_is_valid(x));
	return *this;
}

BigIntClass& BigIntClass::operator=(USmallInt v) {
	bigint_set_usmall(&x, v);
	assert(bigint_is_valid(x));
	return *this;
}

BigIntClass& BigIntClass::operator=(signed v) {
	bigint_set_small(&x, v);
	assert(bigint_is_valid(x));
	return *this;
}

BigIntClass& BigIntClass::operator=(unsigned v) {
	bigint_set_usmall(&x, v);
	assert(bigint_is_valid(x));
	return *this;
}

bool BigIntClass::operator==(const BigIntClass& b) const {
	return bigint_cmp(x, b.x) == 0;
}

bool BigIntClass::operator!=(const BigIntClass& b) const {
	return bigint_cmp(x, b.x) != 0;
}

bool BigIntClass::operator>(const BigIntClass& b) const {
	return bigint_cmp(x, b.x) > 0;
}

bool BigIntClass::operator<(const BigIntClass& b) const {
	return bigint_cmp(x, b.x) < 0;
}

bool BigIntClass::operator>=(const BigIntClass& b) const {
	return bigint_cmp(x, b.x) >= 0;
}

bool BigIntClass::operator<=(const BigIntClass& b) const {
	return bigint_cmp(x, b.x) <= 0;
}

bool BigIntClass::operator==(SmallInt v) const {
	return bigint_cmp_small(x, v) == 0;
}

bool BigIntClass::operator!=(SmallInt v) const {
	return bigint_cmp_small(x, v) != 0;
}

bool BigIntClass::operator>(SmallInt v) const {
	return bigint_cmp_small(x, v) > 0;
}

bool BigIntClass::operator<(SmallInt v) const {
	return bigint_cmp_small(x, v) < 0;
}

bool BigIntClass::operator>=(SmallInt v) const {
	return bigint_cmp_small(x, v) >= 0;
}

bool BigIntClass::operator<=(SmallInt v) const {
	return bigint_cmp_small(x, v) <= 0;
}

BigIntClass BigIntClass::operator-() const {
	BigIntClass z;
	bigint_neg(x, &z.x);
	assert(bigint_is_valid(z.x));
	return z;
}

BigIntClass BigIntClass::operator+(const BigIntClass& b) const {
	BigIntClass c;
	bigint_add(x, b.x, &c.x);
	assert(bigint_is_valid(c.x));
	return c;
}

BigIntClass BigIntClass::operator-(const BigIntClass& b) const {
	BigIntClass c;
	bigint_sub(x, b.x, &c.x);
	assert(bigint_is_valid(c.x));
	return c;
}

BigIntClass BigIntClass::operator*(const BigIntClass& b) const {
	BigIntClass c;
	bigint_mul(x, b.x, &c.x);
	assert(bigint_is_valid(c.x));
	return c;
}

BigIntClass BigIntClass::operator/(const BigIntClass& b) const {
	BigIntClass c;
	bigint_div(x, b.x, &c.x, NULL);
	assert(bigint_is_valid(c.x));
	return c;
}

BigIntClass BigIntClass::operator%(const BigIntClass& b) const {
	BigIntClass c;
	bigint_div(x, b.x, NULL, &c.x);
	assert(bigint_is_valid(c.x));
	return c;
}

std::ostream& operator << (std::ostream& outs, const BigIntClass& z) {
	BigInt_FormatSpec bifs = { 0 };
	bifs.add_spaces = true;
	bifs.add_prefix = (outs.flags() & std::ios_base::showbase);
	char* str;
	if (outs.flags() & std::ios_base::dec) {
		str = bigint_dec_str(z.x, bifs);
	} else {
		str = bigint_hex_str(z.x, bifs);
	}
	outs << str;
	free(str);
	return outs;
}

std::istream& operator >> (std::istream& ins, BigIntClass& z) {
	std::string s;
	ins >> s;
	bigint_scan(s.data(), s.size(), &z.x);
	assert(bigint_is_valid(z.x));
	return ins;
}

BigIntClass& BigIntClass::operator+=(const BigIntClass& b) {
	bigint_add(x, b.x, &x);
	assert(bigint_is_valid(x));
	return *this;
}

BigIntClass& BigIntClass::operator-=(const BigIntClass& b) {
	bigint_sub(x, b.x, &x);
	assert(bigint_is_valid(x));
	return *this;
}

BigIntClass& BigIntClass::operator*=(const BigIntClass& b) {
	bigint_mul(x, b.x, &x);
	assert(bigint_is_valid(x));
	return *this;
}

BigIntClass& BigIntClass::operator/=(const BigIntClass& b) {
	bigint_div(x, b.x, &x, NULL);
	assert(bigint_is_valid(x));
	return *this;
}

BigIntClass& BigIntClass::operator%=(const BigIntClass& b) {
	bigint_div(x, b.x, NULL, &x);
	assert(bigint_is_valid(x));
	return *this;
}

BigIntClass& BigIntClass::add(const BigIntClass& a, const BigIntClass& b) {
	bigint_add(a.x, b.x, &x);
	assert(bigint_is_valid(x));
	return *this;
}

BigIntClass& BigIntClass::sub(const BigIntClass& a, const BigIntClass& b) {
	bigint_sub(a.x, b.x, &x);
	assert(bigint_is_valid(x));
	return *this;
}

BigIntClass& BigIntClass::mul(const BigIntClass& a, const BigIntClass& b) {
	bigint_mul(a.x, b.x, &x);
	assert(bigint_is_valid(x));
	return *this;
}

BigIntClass& BigIntClass::div(const BigIntClass& a, const BigIntClass& b) {
	bigint_div(a.x, b.x, &x, NULL);
	assert(bigint_is_valid(x));
	return *this;
}

BigIntClass& BigIntClass::mod(const BigIntClass& a, const BigIntClass& b) {
	bigint_div(a.x, b.x, NULL, &x);
	assert(bigint_is_valid(x));
	return *this;
}

BigIntClass& BigIntClass::operator++() {
	bigint_add(x, ONE, &x);
	assert(bigint_is_valid(x));
	return *this;
}

BigIntClass  BigIntClass::operator++(int) {
	BigIntClass copy(*this);
	++(*this);
	assert(bigint_is_valid(x));
	return copy;
}

BigIntClass& BigIntClass::operator--() {
	bigint_sub(x, ONE, &x);
	assert(bigint_is_valid(x));
	return *this;
}

BigIntClass  BigIntClass::operator--(int) {
	BigIntClass copy(*this);
	--(*this);
	assert(bigint_is_valid(x));
	return copy;
}

bool BigIntClass::bit(size_t bitpos) const {
	return bigint_bit_get(x, bitpos);
}

BigIntClass& BigIntClass::set_bit(size_t bitpos) {
	bigint_bit_set(&x, bitpos);
	assert(bigint_is_valid(x));
	return *this;
}

BigIntClass& BigIntClass::unset_bit(size_t bitpos) {
	bigint_bit_unset(x, bitpos);
	assert(bigint_is_valid(x));
	return *this;
}

BigIntClass& BigIntClass::set_bit(size_t bitpos, bool val) {
	bigint_bit_set_to(&x, bitpos, val);
	assert(bigint_is_valid(x));
	return *this;
}

BigIntClass& BigIntClass::toggle_bit(size_t bitpos) {
	bigint_bit_toggle(&x, bitpos);
	assert(bigint_is_valid(x));
	return *this;
}

BigIntClass& BigIntClass::randomize(size_t width) {
	bigint_rand(&x, width);
	assert(bigint_is_valid(x));
	return *this;
}

BigIntClass BigIntClass::random(size_t width) {
	BigIntClass z;
	bigint_rand(&z.x, width);
	assert(bigint_is_valid(z.x));
	return z;
}

BigIntClass BigIntClass::abs() {
	BigIntClass z;
	bigint_abs(x, &z.x);
	assert(bigint_is_valid(z.x));
	return z;
}

BigIntClass BigIntClass::neg() {
	BigIntClass z;
	bigint_neg(x, &z.x);
	assert(bigint_is_valid(z.x));
	return z;
}

BigIntClass& BigIntClass::to_abs() {
	bigint_abs(x, &x);
	assert(bigint_is_valid(x));
	return *this;
}

BigIntClass& BigIntClass::to_neg() {
	bigint_neg(x, &x);
	assert(bigint_is_valid(x));
	return *this;
}

BigIntClass BigIntClass::operator<<(size_t shift) const {
	BigIntClass z;
	bigint_lshift(x, shift, &z.x);
	assert(bigint_is_valid(z.x));
	return z;
}

BigIntClass BigIntClass::operator>>(size_t shift) const {
	BigIntClass z;
	bigint_rshift(x, shift, &z.x);
	assert(bigint_is_valid(z.x));
	return z;
}

BigIntClass& BigIntClass::operator<<=(size_t shift) {
	bigint_lshift(x, shift, &x);
	assert(bigint_is_valid(x));
	return *this;
}

BigIntClass& BigIntClass::operator>>=(size_t shift) {
	bigint_rshift(x, shift, &x);
	assert(bigint_is_valid(x));
	return *this;
}

BigIntClass& BigIntClass::copy(Block* src_blocks, size_t size) {
	bigint_copy_blocks(&x, src_blocks, size);
	assert(bigint_is_valid(x));
	return *this;
}

BigIntClass BigIntClass::from_little_endian(const void* src, size_t bitpos, size_t width) {
	BigIntClass z;
	bigint_from_little_endian(src, bitpos, width, &z.x);
	assert(bigint_is_valid(z.x));
	return z;
}

void BigIntClass::to_little_endian(void* dst, size_t bitpos) const {
	bigint_to_little_endian(x, dst, bitpos);
	assert(bigint_is_valid(x));
}

void BigIntClass::warn_bad_recpr(bool y) {
	bigint_warn_bad_recpr(y);
}

void BigIntClass::memstat() {
	bigint_memstat();
}

void BigIntClass::bad_recpr_stat() {
	bigint_bad_recpr_stat();
}

BigIntClass BigIntClass::recpr(size_t precision) const {
	BigIntClass z;
	bigint_recpr(x, precision, &z.x);
	assert(bigint_is_valid(z.x));
	return z;
}

BigIntClass& BigIntClass::div(const BigIntClass& a, const BigIntClass& b, const BigIntClass& recpr, size_t precision) {
	bigint_div_recpr(a.x, b.x, recpr.x, precision, &x, NULL);
	assert(bigint_is_valid(x));
	return *this;
}

BigIntClass& BigIntClass::mod(const BigIntClass& a, const BigIntClass& b, const BigIntClass& recpr, size_t precision) {
	bigint_div_recpr(a.x, b.x, recpr.x, precision, NULL, &x);
	assert(bigint_is_valid(x));
	return *this;
}

BigIntClass BigIntClass::powmod(const BigIntClass& exp, const BigIntClass& mod) const {
	BigIntClass z;
	bigint_powmod(x, exp.x, mod.x, &z.x);
	assert(bigint_is_valid(z.x));
	return z;
}

BigIntClass& BigIntClass::crypt_randomize(size_t width) {
	bigint_crypt_rand(width, &x);
	assert(bigint_is_valid(x));
	return *this;
}

BigIntClass BigIntClass::crypt_random(size_t width) {
	BigIntClass z;
	bigint_crypt_rand(width, &z.x);
	assert(bigint_is_valid(z.x));
	return z;
}

bool BigIntClass::crypt_init() {
	return ::crypt_init();
}

bool BigIntClass::crypt_finish() {
	return ::crypt_finish();
}

