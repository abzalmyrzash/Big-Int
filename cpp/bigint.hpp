#pragma once

#include <stdint.h>
#include <iostream>

class BigIntClass {
	private:
		typedef struct bigint* BigInt;
		BigInt x;
		static bool is_log;

	public:
		constexpr static int BLOCK_WIDTH = 64;
		typedef uint64_t Block;
		typedef int64_t  SmallInt;
		typedef uint64_t USmallInt;

		static void init();
		static bool crypt_init();
		static void finish();
		static bool crypt_finish();
		static void log(bool y);
		static void memstat();
		static void warn_bad_recpr(bool y);
		static void bad_recpr_stat();

		size_t size() const;
		size_t width() const;
		size_t cap() const;
		const Block* data() const;
		
		static size_t width_to_size(size_t width);

		BigIntClass& copy(Block* src_blocks, size_t size);
		static BigIntClass from_little_endian(const void* src, size_t bitpos, size_t width);
		void to_little_endian(void* dst, size_t bitpos) const;

		void reserve(size_t cap, bool keep_data = true);
		void destroy();

		BigIntClass(const BigIntClass&);
		BigIntClass(BigIntClass&&);
		BigIntClass();
		BigIntClass(SmallInt v);
		BigIntClass(USmallInt v);
		BigIntClass(signed v);
		BigIntClass(unsigned v);

		~BigIntClass();

		BigIntClass& operator=(const BigIntClass& b);
		BigIntClass& operator=(BigIntClass&& b);
		BigIntClass& operator=(SmallInt v);
		BigIntClass& operator=(USmallInt v);
		BigIntClass& operator=(signed v);
		BigIntClass& operator=(unsigned v);

		bool operator==(const BigIntClass& b) const;
		bool operator!=(const BigIntClass& b) const;
		bool operator> (const BigIntClass& b) const;
		bool operator< (const BigIntClass& b) const;
		bool operator>=(const BigIntClass& b) const;
		bool operator<=(const BigIntClass& b) const;

		bool operator==(SmallInt v) const;
		bool operator!=(SmallInt v) const;
		bool operator> (SmallInt v) const;
		bool operator< (SmallInt v) const;
		bool operator>=(SmallInt v) const;
		bool operator<=(SmallInt v) const;

		BigIntClass operator-() const;

		BigIntClass operator+(const BigIntClass& b) const;
		BigIntClass operator-(const BigIntClass& b) const;
		BigIntClass operator*(const BigIntClass& b) const;
		BigIntClass operator/(const BigIntClass& b) const;
		BigIntClass operator%(const BigIntClass& b) const;

		BigIntClass& operator+=(const BigIntClass& b);
		BigIntClass& operator-=(const BigIntClass& b);
		BigIntClass& operator*=(const BigIntClass& b);
		BigIntClass& operator/=(const BigIntClass& b);
		BigIntClass& operator%=(const BigIntClass& b);

		BigIntClass& operator++();
		BigIntClass  operator++(int);
		BigIntClass& operator--();
		BigIntClass  operator--(int);

		BigIntClass operator<<(size_t shift) const;
		BigIntClass operator>>(size_t shift) const;

		BigIntClass& operator<<=(size_t shift);
		BigIntClass& operator>>=(size_t shift);

		BigIntClass& add(const BigIntClass& a, const BigIntClass& b);
		BigIntClass& sub(const BigIntClass& a, const BigIntClass& b);
		BigIntClass& mul(const BigIntClass& a, const BigIntClass& b);
		BigIntClass& div(const BigIntClass& a, const BigIntClass& b);
		BigIntClass& mod(const BigIntClass& a, const BigIntClass& b);

		BigIntClass  recpr(size_t precision) const;
		BigIntClass& div(const BigIntClass& a, const BigIntClass& b, const BigIntClass& recpr, size_t precision);
		BigIntClass& mod(const BigIntClass& a, const BigIntClass& b, const BigIntClass& recpr, size_t precision);

		BigIntClass powmod(const BigIntClass& exp, const BigIntClass& mod) const;

		BigIntClass abs();
		BigIntClass neg();

		BigIntClass& to_abs();
		BigIntClass& to_neg();

		bool bit(size_t bitpos) const;
		BigIntClass& set_bit(size_t bitpos);
		BigIntClass& unset_bit(size_t bitpos);
		BigIntClass& set_bit(size_t bitpos, bool val);
		BigIntClass& toggle_bit(size_t bitpos);

		BigIntClass& randomize(size_t width);
		static BigIntClass random(size_t width);

		BigIntClass& crypt_randomize(size_t width);
		static BigIntClass crypt_random(size_t width);

		friend std::ostream& operator << (std::ostream& outs, const BigIntClass& z);
		friend std::istream& operator >> (std::istream& ins, BigIntClass& z);
};

