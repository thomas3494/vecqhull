FC = gfortran
FLAGS = -O3 -march=native -mtune=native -J${HOME}/.local/include
LFLAGS = -lfvqhull -lvqhull -fopenmp -lm
EXAMPLES = $(patsubst %.f90, %, $(wildcard *.f90))

.PHONY: all clean

all: $(EXAMPLES)

%: %.f90
	$(FC) $(FLAGS) $^ -o $@ $(LFLAGS)

clean:
	$(RM) $(EXAMPLES)
