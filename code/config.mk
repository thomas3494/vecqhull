CXX = clang++
CC = clang

# Do not use -Ofast / -ffast-math as this may break the evaluation of orient(p, q, p) to exact zero!
FLAGS = -O3 -march=native -mtune=native -Wall -Wextra -Wno-unused-command-line-argument -Wno-deprecated -Wno-vla-cxx-extension -fopenmp
DBUG_FLAGS = -fsanitize=address -fsanitize=undefined -ggdb
FLAGS += $(DBUG_FLAGS)
LFLAGS = -lm -fopenmp -lhwy
