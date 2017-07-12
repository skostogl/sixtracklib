#CFLAGS= -O1 -fPIC
#CFLAGS= -g -fPIC -Wall -Werror
CFLAGS= -g -fPIC -Wall -Werror -DVECTOR
#CFLAGS= -std=c99 -g -O3 -fPIC -msse -Wall -Werror
#CFLAGS= -qopt-report=5 -std=c99 -g -O3 -fPIC -Wall -Werror
#CC=icc
CC=gcc
#LIBS=-lm
LIBS=


all: lib

lib: block.so


block.o: block.c block.h beam.h track.c track.h value.h
	$(CC) $(CFLAGS) -c -o $@ $<

block_vec.o: block_vec.c block.h beam.h track.c track.h value.h
	$(CC) $(CFLAGS) -c -o $@ $<

block.so: block.o
	$(CC)  -shared $(LIBS) -o $@ $^

clean:
	-rm block.so block.o
