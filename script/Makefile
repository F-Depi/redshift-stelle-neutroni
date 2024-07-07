CC = gcc
CFLAGS = -I${CURDIR}/include

OBJECTS_TEST = test.c fun.c
OBJECTS_MAIN = main.c fun.c

all: test main

test: $(OBJECTS_TEST)
	$(CC) $(CFLAGS) -o test.x $(OBJECTS_TEST) -lm

main: $(OBJECTS_TEST)
	$(CC) $(CFLAGS) -o main.x $(OBJECTS_MAIN) -lm
