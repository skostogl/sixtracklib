DEBFLAGS= -Wall -Werror
OPTFLAGS= -O3
CFLAGS= $(OPTFLAGS) -std=c99 -g -fPIC -Wall -Werror
VECFLAGS= -DVECTOR
#CFLAGS= -std=c99 -g -O3 -fPIC -msse -Wall -Werror
#CFLAGS= -qopt-report=5 -std=c99 -g -O3 -fPIC -Wall -Werror
#CC=icc
CC=gcc
#LIBS=-lm
LIBS=


all: lib

lib: block.so block_soa.so


block.o: block.c block.h beam.h track.c track.h value.h
	$(CC) $(CFLAGS) -c -o $@ $<

block_soa.o: block.c block.h beam.h track.c track.h value.h
	$(CC) $(CFLAGS) $(VECFLAGS) -c -o $@ $<

block.so: block.o
	$(CC)  -shared $(LIBS) -o $@ $^

block_soa.so: block_soa.o
	$(CC)  -shared $(LIBS) -o $@ $^

clean:
	-rm block.so block_soa.so block.o block_soa.o
