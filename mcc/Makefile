CC=gcc
LD=g++
CPP=g++

all:	tmcc

tmcc:	mcc.o
	$(LD) -std=c++17 -O3 mcc.o -o tmcc -fopenmp

mcc.o:	mcc.c
	$(CC) -c mcc.c -O3 -fopenmp

clean:
	rm -f *.o tmcc *~
