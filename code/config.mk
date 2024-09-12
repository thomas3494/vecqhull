# Default, for benchmarking BUILD=RELEASE make
BUILD ?= DEBUG

CXX = g++
CC = gcc

# Do not use -Ofast / -ffast-math as this may break the evaluation of orient(p, q, p) to exact zero!
FLAGS = -O3 -march=native -mtune=native -Wall -Wextra -Wno-unused-command-line-argument -Wno-deprecated -fopenmp
DEBUG_FLAGS = -fsanitize=address -fsanitize=undefined -ggdb
RELEASE_FLAGS = -D NDEBUG
FLAGS += $(${BUILD}_FLAGS)
LFLAGS = -lm -fopenmp -lhwy
