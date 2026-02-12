#include <stdint.h>
#include <iostream>

typedef struct bigint* BigInt;
typedef int64_t SmallInt;
typedef uint64_t USmallInt;

class BigIntClass {
	public:
		BigInt x;
		BigIntClass();
		BigIntClass(SmallInt v);
		~BigIntClass();
		BigIntClass operator=(SmallInt v);
		BigIntClass operator+(const BigIntClass& b);
		BigIntClass operator-(const BigIntClass& b);
		BigIntClass operator*(const BigIntClass& b);
		BigIntClass operator/(const BigIntClass& b);
		BigIntClass operator%(const BigIntClass& b);
};

std::ostream& operator << (std::ostream& outs, const BigIntClass& z);
std::istream& operator >> (std::istream& ins, BigIntClass& z);

