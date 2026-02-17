#include <string>

void initBase64();

void encodeBase64(const std::string& str, std::string& base64);

void decodeBase64(const std::string& base64, std::string& decoded);

void encodeBase64(const void* str, size_t len, std::string& base64);

void decodeBase64(const char* base64, size_t len, std::string& decoded);
