CC = gcc48
CFLAGS= -I/usr/local/include -g -Wunused -Wreturn-type -std=c99
CPPFLAGS= -I/usr/local/include -g -Wunused -Wreturn-type
LDFLAGS = -L/usr/local/lib -lm -lc -lgsl -lgslcblas
CPP = c++

all: psphere.c psphere.h
	$(CC) $(CFLAGS) -c  psphere.c
	$(CC) $(LDFLAGS) -o psphere psphere.o 
