CC = gcc
CXX = g++
FLAGS = -O3 -march=native -mtune=native
LFLAGS = -lvqhull -fopenmp -lm
EXAMPLES = $(filter-out util, $(patsubst %.c, %, $(wildcard *.c)))

.PHONY: all clean

all: $(EXAMPLES)

%: %.c util.o
	$(CC) $(FLAGS) $^ -o $@ $(LFLAGS)

util.o: util.c
	$(CC) $(FLAGS) -c $< -o $@

# For paper, should probably remove when releasing the library
cgal_bench: cgal_bench.cpp
	$(CXX) $(FLAGS) $< -o $@ -lgmp

clean:
	$(RM) $(EXAMPLES) util.o
