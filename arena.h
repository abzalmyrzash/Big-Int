#include "utils.h"

#define ARENA_BASE_POS (sizeof(mem_arena))
#define ARENA_ALIGN (sizeof(void*))

typedef struct {
    size_t reserve_size;
    size_t commit_size;

    size_t pos;
    size_t commit_pos;
} mem_arena;

mem_arena* arena_create(size_t reserve_size, size_t commit_size);
void arena_destroy(mem_arena* arena);
void* arena_push(mem_arena* arena, size_t size, b32 non_zero);
void arena_pop(mem_arena* arena, size_t size);
void arena_pop_to(mem_arena* arena, size_t pos);
void arena_clear(mem_arena* arena);

#define PUSH_STRUCT(arena, T) (T*)arena_push((arena), sizeof(T), false)
#define PUSH_STRUCT_ZERO(arena, T) (T*)arena_push((arena), sizeof(T), true)
#define PUSH_ARRAY(arena, T, n) (T*)arena_push((arena), sizeof(T) * (n), false)
#define PUSH_ARRAY_ZERO(arena, T, n) (T*)arena_push((arena), sizeof(T) * (n), true)
#define POP_ARRAY(arena, T, n) arena_pop((arena), sizeof(T) * (n))

u32 plat_get_pagesize(void);

void* plat_mem_reserve(size_t size);
b32 plat_mem_commit(void* ptr, size_t size);
b32 plat_mem_decommit(void* ptr, size_t size);
b32 plat_mem_release(void* ptr, size_t size);

