#CFLAGS= -O1 -fPIC
#CFLAGS= -g -fPIC -Wall -Werror
#CFLAGS= -std=c99 -g -O3 -fPIC -msse -Wall -Werror
CFLAGS= -qopt-report=5 -std=c99 -g -O3 -fPIC -xAVX2 -qopenmp -Wall -Werror
CC=icc
#LIBS=-lm
LIBS=-qopenmp


all: lib

lib: block.so


block.o: block.c block.h beam.h track.c track.h value.h
	$(CC) $(CFLAGS) -c -o $@ $<

block.so: block.o
	$(CC)  -shared $(LIBS) -o $@ $^

clean:
	-rm block.so block.o
