#pragma once

#define ELOG(...) fprintf(stderr, ##__VA_ARGS__)
#define ELOG_ERROR() ELOG("ERROR in %s:%d:%s\n", __FILE__, __LINE__, __func__)
#define ELOG_STR(STR, ...) ELOG("ERROR in %s:%d:%s: " STR "\n", __FILE__, __LINE__, __func__, ##__VA_ARGS__)
#define ELOG_CODE(CODE) ELOG("ERROR %d in %s:%d:%s: %s\n", CODE, __FILE__, __LINE__, __func__, strerror(CODE))
#define ELOG_CODE_STR(CODE, STR, ...) ELOG("ERROR %d in %s:%d:%s: %s; " STR "\n", CODE, __FILE__, __LINE__, __func__, strerror(CODE), ##__VA_ARGS__)

#define WLOG(...) fprintf(stderr, ##__VA_ARGS__)
#define WLOG_WARN() WLOG("WARNING in %s:%d:%s\n", __FILE__, __LINE__, __func__)
#define WLOG_STR(STR, ...) WLOG("WARNING in %s:%d:%s: " STR "\n", __FILE__, __LINE__, __func__, ##__VA_ARGS__)
