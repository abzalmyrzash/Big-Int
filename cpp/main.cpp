#include "bigint.hpp"
extern "C" {
	#include "..\bigint.h"
};

using namespace std;

int main() {
	bigint_init();
	bigint_log_memory(true);
	BigIntClass a, b;
	cin >> a >> b;
	cout << (a + b) * b << endl;
	cout << b - a << endl;
}
