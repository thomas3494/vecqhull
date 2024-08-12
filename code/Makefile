# TODO build system is sloppy, should handle the C++ / C distinction properly

include config.mk

.PHONY: all clean

all: bin/leftright_test bin/test_partition bin/test_partitionp bin/test_quickhull bin/index2points bin/text2bin bin/test_quickhull_par

bin/common.o: src/common.cpp src/common.h
	$(CXX) $(FLAGS) -c $< -o $@ $(LFLAGS)

bin/quickhull.o: src/quickhull.cpp 
	$(CXX) $(FLAGS) -c $^ -o $@ $(LFLAGS)

bin/leftright_test: src/leftright_test.c bin/common.o
	$(CXX) $(FLAGS) $^ -o $@ $(LFLAGS)

bin/test_partition: src/test_partition.c bin/common.o
	$(CXX) $(FLAGS) $^ -o $@ $(LFLAGS)

bin/test_partitionp: src/test_partitionp.c bin/common.o
	$(CXX) $(FLAGS) $^ -o $@ $(LFLAGS)

bin/test_quickhull: src/test_quickhull.c bin/quickhull.o bin/common.o
	$(CXX) $(FLAGS) $^ -o $@ $(LFLAGS)

bin/test_quickhull_par: src/test_quickhull_par.c bin/quickhull.o bin/common.o
	$(CXX) $(FLAGS) $^ -o $@ $(LFLAGS)

bin/index2points: src/index2points.c
	$(CXX) $(FLAGS) $^ -o $@ $(LFLAGS)

bin/text2bin: src/text2bin.c
	$(CXX) $(FLAGS) $^ -o $@ $(LFLAGS)

clean:
	$(RM) bin/*
