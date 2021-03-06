CPP=g++
SRCS = $(wildcard *.cc)

PROGS = $(patsubst %.cc,%,$(SRCS))

all: $(PROGS)

%: %.cc
	$(CPP) -O3 -o $@ $< -lgomp -fopenmp
clean:
	rm -f $(PROGS)

