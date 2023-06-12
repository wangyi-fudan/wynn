all:	benchmark
benchmark:	benchmark.cpp wynn.hpp makefile
	g++ benchmark.cpp -o benchmark -Ofast -fopenmp -march=native -Wall

