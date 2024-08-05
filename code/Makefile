include config.mk

.PHONY: all clean

all: bin/leftright_test bin/argmax_bench bin/leftright_test_soa

bin/common.o: src/common.cpp src/common.h
	$(CXX) $(FLAGS) -c $< -o $@ $(LFLAGS)

bin/leftright_test: src/leftright_test.c bin/common.o
	$(CXX) $(FLAGS) $^ -o $@ $(LFLAGS)

bin/leftright_test_soa: src/leftright_test_soa.c bin/common.o
	$(CXX) $(FLAGS) $^ -o $@ $(LFLAGS)

bin/argmax_bench: src/argmax_bench.c bin/common.o
	$(CXX) $(FLAGS) $^ -o $@ $(LFLAGS)

clean:
	$(RM) bin/*
