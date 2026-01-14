#include "arena.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

mem_arena* arena_create(size_t reserve_size, size_t commit_size) {
    u32 pagesize = plat_get_pagesize();

    reserve_size = ALIGN_UP_POW2(reserve_size, pagesize);
    commit_size = ALIGN_UP_POW2(commit_size, pagesize);

    mem_arena* arena = plat_mem_reserve(reserve_size);

    if (!plat_mem_commit(arena, commit_size)) {
        return NULL;
    }

    arena->reserve_size = reserve_size;
    arena->commit_size = commit_size;
    arena->pos = ARENA_BASE_POS;
    arena->commit_pos = commit_size;
    
    return arena;
}

void arena_destroy(mem_arena* arena) {
    plat_mem_release(arena, arena->reserve_size);
}

void* arena_push(mem_arena* arena, size_t size, b32 zero) {
    size_t pos_aligned = ALIGN_UP_POW2(arena->pos, ARENA_ALIGN);
    size_t new_pos = pos_aligned + size;

    if (new_pos > arena->reserve_size) { 
		fprintf(stderr, "ERROR: out of memory in arena_push\n");
		abort();
		return NULL;
	}

    if (new_pos > arena->commit_pos) {
        size_t new_commit_pos = new_pos;
        new_commit_pos += arena->commit_size - 1;
        new_commit_pos -= new_commit_pos % arena->commit_size;
        new_commit_pos = MIN(new_commit_pos, arena->reserve_size);

        u8* mem = (u8*)arena + arena->commit_pos;
        size_t commit_size = new_commit_pos - arena->commit_pos;

        if (!plat_mem_commit(mem, commit_size)) {
			fprintf(stderr, "ERROR: out of memory in arena_push\n");
			abort();
            return NULL;
        }

        arena->commit_pos = new_commit_pos;
    }

    arena->pos = new_pos;

    u8* out = (u8*)arena + pos_aligned;

    if (zero) {
        memset(out, 0, size);
    }

    return out;
}

void arena_pop(mem_arena* arena, size_t size) {
    size = MIN(size, arena->pos - ARENA_BASE_POS);
    arena->pos -= size;
}

void arena_pop_to(mem_arena* arena, size_t pos) {
    size_t size = pos < arena->pos ? arena->pos - pos : 0;
    arena_pop(arena, size);
}

void arena_clear(mem_arena* arena) {
    arena_pop_to(arena, ARENA_BASE_POS);
}

#if defined(_WIN32)

#define WIN32_LEAN_AND_MEAN
#include <windows.h>

u32 plat_get_pagesize(void) {
    SYSTEM_INFO sysinfo = { 0 };
    GetSystemInfo(&sysinfo);

    return sysinfo.dwPageSize;
}

void* plat_mem_reserve(size_t size) {
    return VirtualAlloc(NULL, size, MEM_RESERVE, PAGE_READWRITE);
}

b32 plat_mem_commit(void* ptr, size_t size) {
    void* ret = VirtualAlloc(ptr, size, MEM_COMMIT, PAGE_READWRITE);
    return ret != NULL;
}

b32 plat_mem_decommit(void* ptr, size_t size) {
    return VirtualFree(ptr, size, MEM_DECOMMIT);
}

b32 plat_mem_release(void* ptr, size_t size) {
    return VirtualFree(ptr, size, MEM_RELEASE);
}


#elif defined(__linux__)

#define _DEFAULT_SOURCE
#define __USE_MISC

#include <unistd.h>
#include <sys/mman.h>

u32 plat_get_pagesize(void) {
    return (u32)sysconf(_SC_PAGESIZE);
}

void* plat_mem_reserve(size_t size) {
    void* out = mmap(NULL, size, PROT_NONE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    if (out == MAP_FAILED) {
        return NULL;
    }
    return out;
}

b32 plat_mem_commit(void* ptr, size_t size) {
    i32 ret = mprotect(ptr, size, PROT_READ | PROT_WRITE);
    return ret == 0;
}

b32 plat_mem_decommit(void* ptr, size_t size) {
    i32 ret = mprotect(ptr, size, PROT_NONE);
    if (ret != 0) return false;
    ret = madvise(ptr, size, MADV_DONTNEED);
    return ret == 0;
}

b32 plat_mem_release(void* ptr, size_t size) {
    i32 ret = munmap(ptr, size);
    return ret == 0;
}

#endif
