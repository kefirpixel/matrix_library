CC:=gcc
CFLAGS:=-Wall -Wextra -Werror -std=c11
GCOVFLAGS:=-fprofile-arcs -ftest-coverage
LDFLAGS:=-lcheck -lgcov -fprofile-arcs --coverage
ifeq ($(shell uname), Linux)
LDFLAGS +=-lpthread -lrt -lm -lsubunit
OPEN_CMD:=xdg-open
else ifeq ($(shell uname), Darwin)
OPEN_CMD:=open
endif

SRC:= $(wildcard s21_*.c)
OBJ = $(addprefix obj/,$(SRC:.c=.o))
OBJ_GCOV = $(addprefix gcov_obj/,$(SRC:.c=.o))
TEST_SRC:= $(wildcard tests/test_*.c)
TEST_OBJ = $(addprefix tests/,$(TEST_SRC:.c=.o))

all: s21_matrix.a test

s21_matrix.a: $(OBJ)
	ar -src $@ $(OBJ)

s21_matrix_gcov.a: $(OBJ_GCOV)
	ar -src $@ $(OBJ_GCOV) 

test: tests/test_s21_matrix.o s21_matrix_gcov.a
	$(CC) -g $< s21_matrix_gcov.a -o $@ $(LDFLAGS) $(GCOVFLAGS)
	clear
	./test

clean:
	rm -f */*.o gcov_obj/*.gc* ./test *.a
	rm -f -r out gcov_obj/*.info .DS_Store
	rm -rf test.dSYM tests/*.o ../.DS_Store
	rm -rf obj/*

gcov_obj/%.o: %.c
	$(CC) $(CFLAGS) $(GCOVFLAGS) -c $< -o $@

obj/%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

gcov_report: test
	lcov -c -d gcov_obj/. -o gcov_obj/coverage.info
	genhtml gcov_obj/coverage.info --output-directory out
	$(OPEN_CMD) out/index.html

rebuild: clean all

clang:
	@cp ../materials/linters/.clang-format .
	clang-format -i *.c *.h tests/*.c
	@rm -f .clang-format

style:
	@cp ../materials/linters/.clang-format .
	clang-format -n *.c *.h tests/*.c
	@rm -f .clang-format

leaks: rebuild
	leaks -atExit -- ./test

valgrind: rebuild
	CK_FORK=no valgrind --vgdb=no --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose ./test