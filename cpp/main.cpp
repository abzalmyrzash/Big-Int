#include "bigint.hpp"
extern "C" {
	#include "..\bigint.h"
};

using namespace std;

int main() {
	bigint_init();
	bigint_log_memory(true);
	BigIntClass a, b;
	char op;
	cin >> a >> op >> b;
	cout << a << op << b << endl;
	BigIntClass c;
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
			c = a / b;
			break;
		case '%':
			c = a % b;
			break;
		default:
			return -1;
	}
	cout << c << endl;
}

