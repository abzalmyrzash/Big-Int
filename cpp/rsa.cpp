#include "bigint.hpp"
#include "bigint_math.hpp"
#include "rsa.hpp"
#include "..\utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <inttypes.h>

int Jacobi(const BigIntClass& a, const BigIntClass& b) {
	if (a == 1) return 1;
	if (a.bit(0)) {
		BigIntClass num = (a - 1) * (b - 1) / 4;
		int sign = (num.bit(0)) ? -1 : 1;
		return Jacobi(b % a, a) * sign;
	}
	BigIntClass num = (b * b - 1) / 8;
	int sign = (num.bit(0)) ? -1 : 1;
	return Jacobi(a / 2, b) * sign;
}

bool primalityTest(const BigIntClass& b) {
	const int pass = 100;
	int cnt = 0;
	size_t b_width = b.width();
	while (cnt < pass) {
		BigIntClass a = BigIntClass::random(b_width);
		if (a >= b) a.unset_bit(b_width - 1);
		if (gcd(a, b) != 1) break;
		int J = Jacobi(a, b);
		BigIntClass c = a.powmod((b - 1) / 2, b);
		if ((b + J) % b == c) {
			cnt++;
		} else {
			break;
		}
	}
	return cnt == pass;
}

void generateKeys(const size_t width, BigIntClass& n, BigIntClass& d, BigIntClass& e, bool printVal) {
	BigIntClass p, q, phi;
	const size_t p_q_diff = 8;

	int cnt;

	std::cout << "Generating p..." << std::endl;
	cnt = 0;
	do {
		p.crypt_randomize(width / 2 + p_q_diff).set_bit(0);
		cnt++;
	} while(!primalityTest(p));

	if (printVal) {
		std::cout << "p = " << std::hex << p << std::endl << std::endl;
		std::cout << "Attempts: " << std::dec << cnt << std::endl << std::endl;
	}

	std::cout << "Generating q..." << std::endl;
	cnt = 0;
	do {
		q.crypt_randomize(width / 2 - p_q_diff).set_bit(0);
		cnt++;
	} while(p == q || !primalityTest(q));

	if (printVal) {
		std::cout << "q = " << std::hex << q << std::endl << std::endl;
		std::cout << "Attempts: " << std::dec << cnt << std::endl << std::endl;
	}

	BigIntClass& max_p_q = p > q ? p : q;
	n = p * q;
	phi = (p - 1) * (q - 1);

	if (printVal) {
		std::cout << "phi = " << std::hex << phi << std::endl << std::endl;
		std::cout << "n = " << std::hex << n << std::endl << std::endl;
	}

	cnt = 0;
	do {
		d.crypt_randomize(width).set_bit(0);
		cnt++;
	} while (!(d > max_p_q && d < phi && gcd(d, phi) == 1));

	if (printVal) {
		std::cout << "d = " << std::hex << d << std::endl << std::endl;
		std::cout << "cnt: " << std::dec << cnt << std::endl << std::endl;
	}

	BigIntClass x0, x1, x2, y;
	BigIntClass a0, b0, a1, b1, a2, b2;

	x0 = phi;
	a0 = 1;
	b0 = 0;
	x1 = d;
	a1 = 0;
	b1 = 1;

	while (1) {
		y.div(x0, x1);

		x2.sub(x0, y * x1);
		if (x2 == 0) break;

		a2.sub(a0, y * a1);
		b2.sub(b0, y * b1);
		/*if (b2 > (phi / 2) || -b2 > phi / 2) {
			printf("EXCEED\n");
		}*/

		x0 = std::move(x1);
		x1 = std::move(x2);
		a0 = std::move(a1);
		a1 = std::move(a2);
		b0 = std::move(b1);
		b1 = b2;
	}

	if (b2 < 0) b2 += phi;
	e = b2;

	if (printVal) {
		std::cout << "e = " << std::hex << e << std::endl << std::endl;
	}
}

void divide(const std::string& str, size_t bitsPerBlock,
		std::vector<BigIntClass>& out) {
	divide(str.data(), str.size(), bitsPerBlock, out);
}

void divide(const char* str, const size_t strLen, size_t bitsPerBlock,
		std::vector<BigIntClass>& out)
{
	size_t total_width = strLen * CHAR_BIT;
	size_t full_blocks = total_width / bitsPerBlock;
	size_t rem_bits    = total_width % bitsPerBlock;
	out.reserve(full_blocks + (rem_bits > 0));
	size_t bitpos = 0;
	size_t i = 0;
	for (; i < full_blocks; i++) {
		out.push_back(BigIntClass::from_little_endian(str, bitpos, bitsPerBlock));
		bitpos += bitsPerBlock;
	}
	out.push_back(BigIntClass::from_little_endian(str, bitpos, rem_bits));
}

void crypt(const std::vector<BigIntClass>& blocks, BigIntClass n, BigIntClass e,
		std::vector<BigIntClass>& out)
{
	for (const BigIntClass& b : blocks) {
		out.push_back(b.powmod(e, n));
	}
}

void merge(const std::vector<BigIntClass>& blocks, size_t bitsPerBlock,
		std::string& out)
{
	size_t total_width = blocks.size() * bitsPerBlock;
	out.resize(ceil_div(total_width, CHAR_BIT), 0);
	size_t bitpos = 0;
	for (const BigIntClass& b : blocks) {
		b.to_little_endian(out.data(), bitpos);
		bitpos += bitsPerBlock;
	}
}

bool testKeys(BigIntClass n, BigIntClass d, BigIntClass e) {
	const size_t n_width = n.width();
	BigIntClass a = BigIntClass::random(n_width);
	if (a >= n) a.unset_bit(n_width - 1);
	return a.powmod(e, n).powmod(d, n) == a && a.powmod(d, n).powmod(e, n) == a;
}

