#define _CRT_SECURE_NO_WARNINGS
#include "rsa.hpp"
#include "base64.hpp"
#include "..\bits.h"
extern "C" {
#include "..\files.h"
#include "..\prompt.h"
#include "..\clipboard.h"
}

const char* keynamesFilename = "keynames.txt";

typedef enum {
	NONE = -1,
	PRIVATE,
	PUBLIC
} KeyType;

bool loadKey(KeyType keyType);

void loadKeys();

void saveKeynames();

#define KEY_WIDTH 2048
#define KEY_BYTES (KEY_WIDTH / CHAR_BIT)

#define FILENAME_SIZE 256
#define COMMAND_SIZE  512
#define MAX_CMD_LEN   2    // without arguments
#define MESSAGE_SIZE  4096

#define KEYS_DIR  "keys"
#define KEYS_PATH "keys/"

#define HELP_PATH "help/"

#define PUB_EXT  ".pub"
#define PRIV_EXT ".prv"

#define CMD_GENERATE_KEYS "g"
#define CMD_SET_KEYS      "k"
#define CMD_ENCRYPT       "e"
#define CMD_ENCRYPT_MANY  "em"
#define CMD_ENCRYPT_FILE  "ef"
#define CMD_DECRYPT       "d"
#define CMD_DECRYPT_FILE  "df"
#define CMD_HELP          "h"
#define CMD_QUIT          "q"

char message[MESSAGE_SIZE];
char filename[FILENAME_SIZE];

const int keysPathLen = sizeof(KEYS_PATH) - 1;

char keyFilenameFull[FILENAME_SIZE] = KEYS_PATH;
char* keyFilename = keyFilenameFull + keysPathLen;

char privFilenameFull[FILENAME_SIZE] = KEYS_PATH;
char* privFilename = privFilenameFull + keysPathLen;

char pubFilenameFull[FILENAME_SIZE] = KEYS_PATH;
char* pubFilename = pubFilenameFull + keysPathLen;

const int pubExtLen = sizeof(PUB_EXT) - 1;
const int privExtLen = sizeof(PRIV_EXT) - 1;
const int maxExtLen = pubExtLen > privExtLen ?
						pubExtLen : privExtLen;
const int maxNoExtSize = FILENAME_SIZE - maxExtLen - keysPathLen;

BigIntClass privN = 0;
BigIntClass privD = 0;
BigIntClass pubN  = 0;
BigIntClass pubE  = 0;

int main(int argc, char* argv[]) {
	initBase64();
	BigIntClass::init();
	BigIntClass::crypt_init();
	BigIntClass::warn_bad_recpr(false);

	char cmd[COMMAND_SIZE];

	int mkdirRes = mkdir(KEYS_DIR);
	if (mkdirRes == DIR_ALREADY_EXISTS)
	{
		// printf("Loading keys...\n");
		loadKeys();
	}

	else if (mkdirRes == DIR_ERROR) {
		printf("Failed to create keys directory\n");
		return 1;
	}

	printf("Welcome to Abzal's %d-bit RSA encryption program!\n"
			"\n"
			"Type h to get help.\n"
			"\n", KEY_WIDTH);

promptCmd:
	printf("> ");

	char c;
	while ((c = getchar()) == ' ');
	ungetc(c, stdin);

	fgets(cmd, COMMAND_SIZE, stdin);
	const int totalCmdLen = strlen(cmd);
	const int cmdLen = strcspn(cmd, " \n");

	if (cmdLen > MAX_CMD_LEN) {
		printf("Too long.\n");
		goto promptCmd;
	} else if (cmdLen == 0) {
		goto promptCmd;
	}

	cmd[cmdLen] = '\0';
	char* args = cmd + cmdLen + 1;
	const int argsLen = strcspn(args, "\n");
	args[argsLen] = '\0';

	if (strcmp(cmd, CMD_GENERATE_KEYS) == 0)
	{
		BigIntClass n, d, e;
		retryGen:
			printf("Generating keys...\n");
			generateKeys(KEY_WIDTH, n, d, e);
			if (!testKeys(n, d, e)) {
				printf("Test failed. Trying generation again.\n");
				goto retryGen;
			}

		printf("Done!\n");

		if (argsLen > 0 && argsLen < maxNoExtSize) {
			memcpy(keyFilename, args, argsLen + 1);
			goto skipKeyName;
		} else {
			if (argsLen >= maxNoExtSize) {
				printf("Name too long!\n");
			}
		}

	promptKeyName:
		printf("Enter key name: ");
		fgets(keyFilename, maxNoExtSize, stdin);

	skipKeyName:
		const int noExtLen = strcspn(keyFilename, ".\n");
		if (noExtLen == 0) {
			goto promptKeyName;
		}

		keyFilename[noExtLen] = '\0';
		strcat(keyFilename, PUB_EXT);

		FILE* file;
		file = fopen(keyFilenameFull, "r");

		if (file) {
			fclose(file);
			printf("%s exists! ", keyFilename);
			printf("Confirm overwrite? (Y/N) ");
			char c = getOption();
			if (c != 'Y') {
				goto promptKeyName;
			}
		}

		file = fopen(keyFilenameFull, "wb");

		if (!file) {
			printf("Could not open %s. Please retry.\n",
					keyFilename);
			goto promptKeyName;
		}

		std::string le1(KEY_BYTES, 0);
		std::string le2(KEY_BYTES, 0);
		n.to_little_endian(le1.data(), 0);
		e.to_little_endian(le2.data(), 0);
		fwrite(le1.data(), KEY_BYTES, 1, file);
		fwrite(le2.data(), KEY_BYTES, 1, file);
		fclose(file);
		le2.assign(KEY_BYTES, 0);

		printf("Public key saved in %s\n", keyFilename);

		keyFilename[noExtLen] = '\0';
		strcat(keyFilename, PRIV_EXT);

		file = fopen(keyFilenameFull, "r");

		if (file) {
			fclose(file);
			printf("%s exists! ", keyFilename);
			printf("Confirm overwrite? (Y/N) ");
			char c = getOption();
			if (c != 'Y') {
				goto promptCmd;
			}
		}

		file = fopen(keyFilenameFull, "wb");

		if (!file) {
			printf("Could not open %s. Please retry.\n",
					keyFilename);
			goto promptCmd;
		}

		d.to_little_endian(le2.data(), 0);
		fwrite(le1.data(), KEY_BYTES, 1, file);
		fwrite(le2.data(), KEY_BYTES, 1, file);
		fclose(file);

		printf("Private key saved in %s\n", keyFilename);

		printf("Set %s as your default private key? (Y/N) ",
				keyFilename);
		char c = getOption();
		if (c == 'Y' || c == 'y') {
			privN = n;
			privD = d;
			const int totalLen = noExtLen + privExtLen;
			memcpy(privFilename, keyFilename, totalLen + 1);
			saveKeynames();
			printf("Saved.\n");
		}

		printf("Remember: only share your public key file, "
				"NOT your private key!\n");
	}

	else if (strcmp(cmd, CMD_SET_KEYS) == 0)
	{
		char* filename = strtok(args, " \n");
		if (filename == NULL) {
			printf("Private key: %s\n", privFilename);
			printf("Public key:  %s\n", pubFilename);
			goto promptCmd;
		}

		KeyType keyType = NONE;
		int cnt = 0;
		while (filename != NULL)
		{
			const char* ext = getExtension(filename);

			if (ext == NULL) {
				printf("No extension!\n");
				goto promptCmd;
			}

			if (strcmp(ext, PUB_EXT) == 0) {
				if (keyType != PUBLIC) {
					keyType = PUBLIC;
				} else {
					printf("Select one public and "
							"one private key!\n");
					goto promptCmd;
				}
			}

			else if (strcmp(ext, PRIV_EXT) == 0) {
				if (keyType != PRIVATE) {
					keyType = PRIVATE;
				} else {
					printf("Select one public and "
							"one private key!\n");
					goto promptCmd;
				}
			}

			else {
				printf("Wrong extension!\n");
				goto promptCmd;
			}

			strcpy(keyFilename, filename);

			loadKey(keyType);

			if (++cnt == 2) break;
			filename = strtok(NULL, " \n");
		}
	}

	else if (strcmp(cmd, CMD_ENCRYPT) == 0 ||
			 strcmp(cmd, CMD_ENCRYPT_MANY) == 0 ||
			 strcmp(cmd, CMD_ENCRYPT_FILE) == 0)
	{
		if (pubN == 0 || pubE == 0) {
			printf("Public key not set!\n");
			goto promptCmd;
		}

		char* str;
		size_t len;
		bool isFile;

		if (strcmp(cmd, CMD_ENCRYPT) == 0) {
			isFile = false;
			if (argsLen == 0) {
				printf("Please provide text!\n");
				goto promptCmd;
			}
			str = args;
			len = argsLen;
		}

		else if (strcmp(cmd, CMD_ENCRYPT_MANY) == 0) {
			isFile = false;
			char ch;
			int cnt;

			if (argsLen > 0) {
				memcpy(message, args, argsLen + 1);
				message[argsLen] = '\n';
				cnt = argsLen + 1;
			} else cnt = 0;

			while (cnt < MESSAGE_SIZE - 1) {
				ch = getchar();
				if (ch == EOF) {
					if (message[cnt - 1] == '\n') {
						message[--cnt] = '\0';
					}
					break;
				}
				message[cnt++] = ch;
			}
			str = message;
			len = cnt;
		}

		else {
			isFile = true;
			FILE* file;

			if (argsLen == 0) {
				printf("Please provide filename!\n");
				goto promptCmd;
			}
			memcpy(filename, args, argsLen + 1);

			file = fopen(filename, "rb");
			if (!file) {
				printf("Could not open file. Please retry.\n");
				goto promptCmd;
			}

			const char* basename = getBasename(filename);
			const int basenameSize = strlen(basename) + 1;

			fseek(file, 0, SEEK_END);
			long fsize = ftell(file);
			fseek(file, 0, SEEK_SET);

			if (fsize == 0) {
				printf("File is empty!\n");
				goto promptCmd;
			}

			len = fsize + basenameSize + sizeof(uint64_t);

			str = (char*)malloc(len);
			encodeLittleEndian64(fsize, str);
			
			fread(str + sizeof(uint64_t), len, 1, file);
			fclose(file);

			memcpy(str + sizeof(uint64_t) + fsize,
					basename, basenameSize);
		}

		printf("Using %s.\n", pubFilename);

		const int bitsPerBlock = pubN.width();
		size_t size;
		std::vector<BigIntClass> blocks;
		std::vector<BigIntClass> encrypted;
		std::string merged;
		std::string base64;

		printf("Dividing into blocks...\n");
		divide(str, len, bitsPerBlock - 1, blocks);
		if (isFile) free(str);

		printf("Encrypting...\n");
		crypt(blocks, pubN, pubE, encrypted);
		blocks.clear();

		printf("Merging...\n");
		merge(encrypted, bitsPerBlock, merged);
		encrypted.clear();

		if (isFile) {
			FILE* file;
		retryEncryptSave:
			printf("Save as: ");
			fgets(filename, FILENAME_SIZE, stdin);
			filename[strcspn(filename, "\n")] = '\0';

			file = fopen(filename, "r");
			if (file) {
				fclose(file);
				printf("%s exists! ", filename);
				printf("Confirm overwrite? (Y/N) ");
				char c = getOption();
				if (c != 'Y') {
					goto retryEncryptSave;
				}
			}

			file = fopen(filename, "wb");
			if (!file) {
				printf("Could not open %s. Please retry.\n",
						filename);
				goto retryEncryptSave;
			}
			
			fwrite(merged.data(), merged.size(), 1, file);
			fclose(file);
			printf("File saved.\n");

			merged.clear();
		}

		else {
			printf("Encoding to Base64...\n");
			encodeBase64(merged, base64);
			merged.clear();
			std::cout << base64 << std::endl;
			copyToClipboard(base64.data());
			printf("Result copied to clipboard automatically.\n");
			base64.clear();
		}
	}

	else if (strcmp(cmd, CMD_DECRYPT) == 0 ||
			 strcmp(cmd, CMD_DECRYPT_FILE) == 0)
	{
		if (privN == 0 || privD == 0) {
			printf("Private key not set!\n");
			goto promptCmd;
		}

		std::string str;
		bool isFile;

		if (strcmp(cmd, CMD_DECRYPT) == 0) {
			isFile = false;
			if (argsLen == 0) {
				printf("Please provide ciphertext!\n");
				goto promptCmd;
			}
			printf("Decoding from Base64...\n");
			decodeBase64(args, argsLen, str);
		}

		else {
			isFile = true;
			FILE* file;

			if (argsLen == 0) {
				printf("Please provide filename!\n");
				goto promptCmd;
			}
			memcpy(filename, args, argsLen + 1);

			file = fopen(filename, "rb");
			if (!file) {
				printf("Could not open file. Please retry.\n");
				goto promptCmd;
			}

			fseek(file, 0, SEEK_END);
			long fsize = ftell(file);
			fseek(file, 0, SEEK_SET);

			if (fsize == 0) {
				printf("File is empty!\n");
				goto promptCmd;
			}

			str.resize(fsize);
			fread(str.data(), fsize, 1, file);
			fclose(file);
		}

		printf("Using %s.\n", privFilename);

		const int bitsPerBlock = privN.width();
		size_t size;
		std::vector<BigIntClass> blocks;
		std::vector<BigIntClass> decrypted;
		std::string merged;

		printf("Dividing into blocks...\n");
		divide(str, bitsPerBlock, blocks);
		if (isFile) str.clear();

		printf("Decrypting...\n");
		crypt(blocks, privN, privD, decrypted);
		blocks.clear();

		printf("Merging...\n");
		merge(decrypted, bitsPerBlock - 1, merged);
		decrypted.clear();

		if (isFile) {
			FILE* file;
			size_t fsize = decodeLittleEndian64(merged.data());

			char* basename = merged.data() + sizeof(uint64_t) + fsize;
			
			printf("Decrypted filename: %s\n", basename);
			printf("Save as %s? (Y/N) ", basename);
			char opt = getOption();
			if (opt == 'Y') {
				strcpy(filename, basename);
				goto skipDecryptSavePrompt;
			}

		decryptSavePrompt:
			printf("Save as: ");
			fgets(filename, FILENAME_SIZE, stdin);
			filename[strcspn(filename, "\n")] = '\0';
		skipDecryptSavePrompt:

			file = fopen(filename, "r");
			if (file) {
				fclose(file);
				printf("%s exists! ", filename);
				printf("Confirm overwrite? (Y/N) ");
				char c = getOption();
				if (c != 'Y') {
					goto decryptSavePrompt;
				}
			}

			file = fopen(filename, "wb");
			if (!file) {
				printf("Could not open %s. Please retry.\n",
						filename);
				goto decryptSavePrompt;
			}

			fwrite(merged.data() + sizeof(size_t), fsize, 1, file);
			fclose(file);
			printf("File saved.\n");

			merged.clear();
		}

		else {
			printf("%s\n", merged.data());
			copyToClipboard(merged.data());
			printf("Result copied to clipboard automatically.\n");
			merged.clear();
		}
	}

	else if (strcmp(cmd, CMD_HELP) == 0)
	{
		int pageNum = 0;

		if (argsLen > 0 && 1 == sscanf(args, "%d", &pageNum)) {
			if (pageNum == 2) goto helpPage2;
			if (pageNum == 3) goto helpPage3;
		}

		if (!printFile(HELP_PATH "1.txt")) {
			printf("Help not found!\n");
			goto promptCmd;
		}

		if (pageNum == 1) goto promptCmd;
		pressAnyKeyToContinue();

	helpPage2:
		if (!printFile(HELP_PATH "2.txt")) {
			printf("Help not found!\n");
			goto promptCmd;
		}

		if (pageNum == 2) goto promptCmd;
		pressAnyKeyToContinue();

	helpPage3:
		if (!printFile(HELP_PATH "3.txt")) {
			printf("Help not found!\n");
			goto promptCmd;
		}
	}

	else if (strcmp(cmd, CMD_QUIT) == 0)
	{
		BigIntClass::finish();
		BigIntClass::crypt_finish();
		return 0;
	}

	else {
		printf("Invalid command!\n");
		goto promptCmd;
	}

	goto promptCmd;
}

void saveKeynames() {
	FILE* file = fopen(keynamesFilename, "w");
	fprintf(file, "%s\n%s", privFilename, pubFilename);
	fclose(file);
}

bool loadKey(KeyType keyType) {
	FILE* file = fopen(keyFilenameFull, "rb");

	if (!file) {
		printf("Could not open file.\n");
		return false;
	}

	fseek(file, 0, SEEK_END);
	long fsize = ftell(file);
	fseek(file, 0, SEEK_SET);
	if (fsize != KEY_BYTES * 2) {
		printf("Failed to load key file. Size invalid.\n");
		return false;
	}

	BigIntClass key1, key2;

	std::string str;
	str.resize(KEY_BYTES);
	fread(str.data(), KEY_BYTES, 1, file);
	key1 = BigIntClass::from_little_endian(str.data(), 0, KEY_WIDTH);
	fread(str.data(), KEY_BYTES, 1, file);
	key2 = BigIntClass::from_little_endian(str.data(), 0, KEY_WIDTH);
	fclose(file);

	if (keyType == PUBLIC) {
		pubN = std::move(key1);
		pubE = std::move(key2);
		strcpy(pubFilename, keyFilename);
		printf("%s successfully loaded.\n", keyFilename);
		saveKeynames();
		return true;
	}
	else {
		privN = std::move(key1);
		privD = std::move(key2);
		strcpy(privFilename, keyFilename);
		printf("%s successfully loaded.\n", keyFilename);
		saveKeynames();
		return true;
	}

}

void loadKeys() {
	FILE* kFile = fopen(keynamesFilename, "r");
	if (!kFile) {
		// printf("Failed to open keynames file.\n");
		return;
	}

	// private key

	if (!fgets(privFilename, FILENAME_SIZE, kFile)) return;
	privFilename[strcspn(privFilename, "\n")] = '\0';
	// printf("Loading %s\n", privFilename);

	std::string str;
	str.resize(KEY_BYTES);
	FILE* file;
	long fsize;

	file = fopen(privFilenameFull, "rb");
	if (!file) {
		printf("Failed to open key file %s.", privFilename);
		return;
	}

	fseek(file, 0, SEEK_END);
	fsize = ftell(file);
	fseek(file, 0, SEEK_SET);
	if (fsize != KEY_BYTES * 2) {
		printf("Failed to load key file %s. Size invalid.\n", privFilename);
		return;
	}

	fread(str.data(), KEY_BYTES, 1, file);
	privN = BigIntClass::from_little_endian(str.data(), 0, KEY_WIDTH);
	fread(str.data(), KEY_BYTES, 1, file);
	privD = BigIntClass::from_little_endian(str.data(), 0, KEY_WIDTH);
	fclose(file);

	// public key

	if (!fgets(pubFilename, FILENAME_SIZE, kFile)) return;
	pubFilename[strcspn(pubFilename, "\n")] = '\0';
	// printf("Loading %s\n", pubFilename);

	fclose(kFile);

	file = fopen(pubFilenameFull, "rb");
	if (!file) {
		printf("Failed to open key file %s.", pubFilename);
		return;
	}

	fseek(file, 0, SEEK_END);
	fsize = ftell(file);
	fseek(file, 0, SEEK_SET);
	if (fsize != KEY_BYTES * 2) {
		printf("Failed to load key file %s. Size invalid.\n", pubFilename);
		return;
	}

	fread(str.data(), KEY_BYTES, 1, file);
	pubN = BigIntClass::from_little_endian(str.data(), 0, KEY_WIDTH);
	fread(str.data(), KEY_BYTES, 1, file);
	pubE = BigIntClass::from_little_endian(str.data(), 0, KEY_WIDTH);
	fclose(file);
}

