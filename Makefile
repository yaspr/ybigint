CC=gcc

#Debug
DFLAGS=-g3

#Shared library
CFLAGS=-fPIC -shared

#Compiler optimzations
OFLAGS=-Ofast -finline-functions -funroll-loops -ftree-vectorize

#Library flags for test
LFLAGS=-lybigint -L.

LIB_IN=bigint.c
LIB_OUT=libybigint.so

TEST_IN=test.c
TEST_OUT=libtest

all: lib test

lib:
	$(CC) $(DFLAGS) $(CFLAGS) $(OFLAGS) $(LIB_IN) -o $(LIB_OUT)

test:
	$(CC) $(DFLAGS) $(LFLAGS) $(OFLAGS) $(TEST_IN) -o $(TEST_OUT)

clean:
	rm -Rf *~ $(LIB_OUT) $(TEST_OUT)
