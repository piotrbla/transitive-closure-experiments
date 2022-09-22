CC=g++
LD=g++
CPP=g++

all:	tlds

tlds:	lds.o
	$(LD) -std=c++17 -O3 lds.o -o tlds -fopenmp

lds.o:	lds.cc
	$(CC) -c lds.cc -O3 -fopenmp

clean:
	rm -f *.o tlds *~
