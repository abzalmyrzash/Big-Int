#include "bigint.hpp"
#include "bigint_math.hpp"
#include "rsa.hpp"
extern "C" {
	#include "..\bigint.h"
};

#include <iostream>
#include <string>
#include <chrono>
#include <cmath>

using namespace std;

BigIntClass fib(unsigned int n) {
	if (n == 0) {
		return 0;
	} else if (n == 1 || n == 2) {
		return 1;
	}

	constexpr double LOG2_GOLDEN_RATIO = 0.69424191363; 
	const size_t cap = ceil(n * LOG2_GOLDEN_RATIO / BIGINT_BLOCK_WIDTH);

	BigIntClass a, b, tmp;
	a.reserve(cap);
	b.reserve(cap);
	tmp.reserve(cap);
	a = 1, b = 1;
	for (int i = 3; i <= n; ++i) {
		tmp = b;
		b += a;
		a = tmp;
	}
	return b;
}

int main() {
	BigIntClass::crypt_init();
	BigIntClass::init();
	BigIntClass::warn_bad_recpr(false);
	while (1) {
		BigIntClass n, d, e;
		cout << "Start" << endl;
		// srand(time(NULL));
		auto start = chrono::high_resolution_clock::now();
		generateKeys(2048, n, d, e);
		auto time1 = chrono::duration<double>(chrono::high_resolution_clock::now() - start);
		cout << std::hex;
		cout << n << endl << endl;
		cout << d << endl << endl;
		cout << e << endl << endl;
		cout << "Time: " <<  time1.count() << endl;
		int bitsPerBlock = n.width() - 1;
		vector<BigIntClass> message;
		message.push_back(BigIntClass::random(bitsPerBlock));
		vector<BigIntClass> encrypted;
		vector<BigIntClass> decrypted;
		crypt(message, n, d, encrypted);
		crypt(encrypted, n, e, decrypted);
		cout << message.back() << endl << endl;
		cout << encrypted.back() << endl << endl;
		cout << decrypted.back() << endl << endl;
		BigIntClass::bad_recpr_stat();
	}
	BigIntClass::finish();
	BigIntClass::crypt_finish();
	return 0;

	/*
	unsigned int fib_n;
	cin >> fib_n;

	auto start = chrono::high_resolution_clock::now();
	BigIntClass ans = fib(fib_n);
	auto time1 = chrono::duration<double>(chrono::high_resolution_clock::now() - start);
	start = chrono::high_resolution_clock::now();
	cout << "Answer: " << std::showbase << ans << endl;
	auto time2 = chrono::duration<double>(chrono::high_resolution_clock::now() - start);
	cout << "Time for fib: " <<  time1.count() << endl;
	cout << "Time to print: " << time2.count() << endl;
	ans.destroy();
	BigIntClass::finish();
	bigint_memstat();
	return 0;
	*/

	BigIntClass a, b, c;
	bool quit;
	while (!quit) {
		char op;
		cin >> a >> op >> b;

		if (a > b) {
			cout << "a > b" << endl;
		} else if (a < b) {
			cout << "a < b" << endl;
		} else {
			cout << "a == b" << endl;
		}
		if (a == 0) int a = 0;

		switch(op) {
			case '+':
				c = a + b;
				break;
			case '-':
				c = a - b;
				break;
			case '*':
				c = a * b;
				break;
			case '/':
				if (b == 0) {
					cout << "CAN'T DIVIDE BY ZERO!" << endl;
					continue;
				}
				c = a / b;
				break;
			case '%':
				if (b == 0) {
					cout << "CAN'T DIVIDE BY ZERO!" << endl;
					continue;
				}
				c = a % b;
				break;
			default:
				bool quit = true;
				break;
		}
		cout << c << endl;
	}
	bigint_memstat();
	bigint_finish();
	return 0;
}

