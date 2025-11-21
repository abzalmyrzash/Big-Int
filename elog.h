#pragma once

#define ELOG(...) fprintf(stderr, __VA_ARGS__)
#define ELOG_ERROR() ELOG("ERROR in %s:%d:%s\n", __FILE__, __LINE__, __FUNCTION__)
#define ELOG_STR(STR) ELOG("ERROR in %s:%d:%s: %s\n", __FILE__, __LINE__, __FUNCTION__, STR)
#define ELOG_CODE(CODE) ELOG("ERROR %d in %s:%d:%s: %s\n", CODE, __FILE__, __LINE__, __FUNCTION__, strerror(CODE))
#define ELOG_CODE_STR(CODE, STR) ELOG("ERROR %d in %s:%d:%s: %s; %s\n", CODE, __FILE__, __LINE__, __FUNCTION__, strerror(CODE), STR)
