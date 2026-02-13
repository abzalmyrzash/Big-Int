#include <stdint.h>
#include <iostream>

typedef struct bigint* BigInt;
typedef int64_t SmallInt;
typedef uint64_t USmallInt;

class BigIntClass {
	private:
		BigInt x;
		int id;
		static bool is_log;

	public:
		static void init();
		static void finish();
		static void log(bool y);

		void reserve(size_t cap, bool keep_data = true);
		void destroy();

		BigIntClass(const BigIntClass&);
		BigIntClass(BigIntClass&&);
		BigIntClass();
		BigIntClass(SmallInt v);

		~BigIntClass();

		BigIntClass& operator=(const BigIntClass& b);
		BigIntClass& operator=(BigIntClass&& b);
		BigIntClass& operator=(SmallInt v);
		BigIntClass& operator=(USmallInt v);
		BigIntClass& operator=(signed v);
		BigIntClass& operator=(unsigned v);

		bool operator==(const BigIntClass& b);
		bool operator>(const BigIntClass& b);
		bool operator<(const BigIntClass& b);
		bool operator>=(const BigIntClass& b);
		bool operator<=(const BigIntClass& b);

		bool operator==(SmallInt v);
		bool operator>(SmallInt v);
		bool operator<(SmallInt v);
		bool operator>=(SmallInt v);
		bool operator<=(SmallInt v);

		BigIntClass operator+(const BigIntClass& b);
		BigIntClass operator-(const BigIntClass& b);
		BigIntClass operator*(const BigIntClass& b);
		BigIntClass operator/(const BigIntClass& b);
		BigIntClass operator%(const BigIntClass& b);

		BigIntClass& operator+=(const BigIntClass& b);
		BigIntClass& operator-=(const BigIntClass& b);
		BigIntClass& operator*=(const BigIntClass& b);
		BigIntClass& operator/=(const BigIntClass& b);
		BigIntClass& operator%=(const BigIntClass& b);

		BigIntClass& operator++();
		BigIntClass  operator++(int);
		BigIntClass& operator--();
		BigIntClass  operator--(int);

		BigIntClass& add(const BigIntClass& a, const BigIntClass& b);
		BigIntClass& sub(const BigIntClass& a, const BigIntClass& b);
		BigIntClass& mul(const BigIntClass& a, const BigIntClass& b);
		BigIntClass& div(const BigIntClass& a, const BigIntClass& b);
		BigIntClass& mod(const BigIntClass& a, const BigIntClass& b);

		friend std::ostream& operator << (std::ostream& outs, const BigIntClass& z);
		friend std::istream& operator >> (std::istream& ins, BigIntClass& z);
};

