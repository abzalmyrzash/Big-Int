CC_FLAGS := -g -O3 --std=c23 -D_CRT_SECURE_NO_WARNINGS

ifdef OS
	RM = del /Q
	FixPath = $(subst /,\,$1)
	CC = clang
	DEBUG := debug
	RELEASE := release
	EXE := main.exe

else
	ifeq ($(shell uname), Linux)
		RM = rm -f
		FixPath = $1
		CC := gcc-14
		LD_FLAGS := -lm -lc
		# CC_FLAGS += -fsanitize=address
		DEBUG := debug-linux
		RELEASE := release-linux
		EXE := main.out
	endif
endif

DEBUG_FLAGS := -DDEBUG
RELEASE_FLAGS := -DNDEBUG

BUILD := $(DEBUG)
EXTRA_FLAGS := $(DEBUG_FLAGS)

CC_FLAGS += $(EXTRA_FLAGS)
CC_FLAGS += $(USER_FLAGS)

.PHONY: all clean debug release

all: $(BUILD)/$(EXE)

debug:
	@echo Debug build
	@$(MAKE) --no-print-directory BUILD=$(DEBUG) EXTRA_FLAGS="$(DEBUG_FLAGS)"

release:
	@echo Release build
	@$(MAKE) --no-print-directory BUILD=$(RELEASE) EXTRA_FLAGS="$(RELEASE_FLAGS)"

$(BUILD)/$(EXE): main.c $(BUILD)/fib.o $(BUILD)/factorial.o $(BUILD)/bigint.o $(BUILD)/bigint_rand.o $(BUILD)/arena.o utils.h
	$(CC) $(CC_FLAGS) main.c $(BUILD)/*.o -o $(BUILD)/$(EXE) $(LD_FLAGS)

$(BUILD)/bigint.o: bigint.c bigint.h bigint_impl.h bigint_impl_basic.h bigint_params.h utils.h
	$(CC) $(CC_FLAGS) -c bigint.c -o $(BUILD)/bigint.o

$(BUILD)/fib.o: fib.c fib.h utils.h
	$(CC) $(CC_FLAGS) -c fib.c -o $(BUILD)/fib.o

$(BUILD)/factorial.o: factorial.c factorial.h utils.h
	$(CC) $(CC_FLAGS) -c factorial.c -o $(BUILD)/factorial.o

$(BUILD)/bigint_rand.o: bigint_rand.c bigint_rand.h utils.h
	$(CC) $(CC_FLAGS) -c bigint_rand.c -o $(BUILD)/bigint_rand.o

$(BUILD)/arena.o: arena.c arena.h utils.h
	$(CC) $(CC_FLAGS) -c arena.c -o $(BUILD)/arena.o

clean:
	$(RM) $(call FixPath,$(DEBUG)/*)
	$(RM) $(call FixPath,$(RELEASE)/*)

