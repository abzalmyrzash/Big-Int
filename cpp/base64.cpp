#include "base64.hpp"
#include <string>
#include "..\utils.h"

static char toBase64[65] =
"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

static char fromBase64[128];

void initBase64() {
	for (int i = 0; i < 64; i++) {
		fromBase64[toBase64[i]] = i;
	}
}

void encodeBase64(const std::string& str, std::string& base64) {
	encodeBase64(str.data(), str.size(), base64);
}

void decodeBase64(const std::string& base64, std::string& decoded) {
	decodeBase64(base64.data(), base64.size(), decoded);
}

void encodeBase64(const void* str, size_t len, std::string& base64) {
	uint8_t* bytes = (uint8_t*)str;
	base64.resize(ceil_div(len * 4, 3), 0);
	size_t i, cnt = 0;

	for (i = 0; i + 2 < len; i += 3) {
		base64[cnt++] = toBase64[bytes[i] >> 2];
		base64[cnt++] = toBase64[((bytes[i] & 3) << 4) |
			(bytes[i + 1] >> 4)];
		base64[cnt++] = toBase64[((bytes[i + 1] & 15) << 2) |
			(bytes[i + 2] >> 6)];
		base64[cnt++] = toBase64[bytes[i + 2] & 63];
	}

	if (i < len) {
		base64[cnt++] = toBase64[bytes[i] >> 2];
		if (i + 1 < len) {
			base64[cnt++] = toBase64[((bytes[i] & 3) << 4) |
				(bytes[i + 1] >> 4)];
			base64[cnt++] = toBase64[(bytes[i + 1] & 15) << 2];
		} else {
			base64[cnt++] = toBase64[(bytes[i] & 3) << 4];
		}
	}

	base64.resize(cnt);
}

void decodeBase64(const char* base64, size_t base64Len, std::string& decoded) {
	decoded.resize(ceil_div(base64Len * 3, 4), 0);
	uint8_t* bytes = (uint8_t*)decoded.data();
	size_t i, cnt = 0;
	uint8_t bix1, bix2, bix3, bix4;

	for (i = 0; i + 3 < base64Len; i += 4) {
		bix1 = fromBase64[base64[i]];
		bix2 = fromBase64[base64[i + 1]];
		bix3 = fromBase64[base64[i + 2]];
		bix4 = fromBase64[base64[i + 3]];
		bytes[cnt++] = (bix1 << 2) | (bix2 >> 4);
		bytes[cnt++] = ((bix2 & 15) << 4) | (bix3 >> 2);
		bytes[cnt++] = ((bix3 & 3) << 6) | bix4;
	}

	if (i < base64Len) {
		bix1 = fromBase64[base64[i]];
		bytes[cnt] = (bix1 << 2);
		if (i + 1 < base64Len) {
			bix2 = fromBase64[base64[i + 1]];
			bytes[cnt++] |= (bix2 >> 4);
			bytes[cnt] = ((bix2 & 15) << 4);
			if (i + 2 < base64Len) {
				bix3 = fromBase64[base64[i + 2]];
				bytes[cnt++] |= (bix3 >> 2);
				bytes[cnt++] = (bix3 & 3 << 6);
			} else cnt++;
		} else cnt++;
	}

	decoded.resize(cnt);
}

