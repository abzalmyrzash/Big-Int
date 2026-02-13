#include "bigint.hpp"
extern "C" {
	#include "..\bigint.h"
};

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
	BigIntClass::init();
	// BigIntClass::log(true);
	unsigned int n;
	cin >> n;

	auto start = chrono::high_resolution_clock::now();
	BigIntClass ans = fib(n);
	auto time1 = chrono::duration<double>(chrono::high_resolution_clock::now() - start);
	start = chrono::high_resolution_clock::now();
	cout << "Answer: " << ans << endl;
	auto time2 = chrono::duration<double>(chrono::high_resolution_clock::now() - start);
	cout << "Time for fib: " <<  time1.count() << endl;
	cout << "Time to print: " << time2.count() << endl;
	ans.destroy();
	BigIntClass::finish();
	bigint_memstat();
	return 0;

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

