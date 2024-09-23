CC = gcc
FLAGS = -O3 -march=native -mtune=native
LFLAGS = -lvqhull -fopenmp -lm -lstdc++
EXAMPLES = $(filter-out util, $(patsubst %.c, %, $(wildcard *.c)))

.PHONY: all clean

all: $(EXAMPLES)

%: %.c util.o
	$(CC) $(FLAGS) $^ -o $@ $(LFLAGS)

util.o: util.c
	$(CC) $(FLAGS) -c $< -o $@

clean:
	$(RM) $(EXAMPLES) util.o
