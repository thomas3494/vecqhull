include config.mk

.PHONY: all clean

all: bin/leftright_test bin/test_partition

bin/common.o: src/common.cpp src/common.h
	$(CXX) $(FLAGS) -c $< -o $@ $(LFLAGS)

bin/leftright_test: src/leftright_test.c bin/common.o
	$(CXX) $(FLAGS) $^ -o $@ $(LFLAGS)

bin/test_partition: src/test_partition.c bin/common.o
	$(CXX) $(FLAGS) $^ -o $@ $(LFLAGS)

clean:
	$(RM) bin/*
