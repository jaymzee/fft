CC = gcc
CXX = g++
CFLAGS = -std=c99 -pedantic-errors -Wall -O3
CXXFLAGS = -std=c++11 -pedantic-errors -Wall -O3
LFLAGS = -lm
COBJS = test_fft benchmark_fft
CXXOBJS = TestFFT TestFFTVector BenchmarkFFT

all: $(COBJS) $(CXXOBJS)

test_fft: test_fft.o fft.o util.o
	$(CC) -o $@ $^ $(LFLAGS)

benchmark_fft: benchmark_fft.o fft.o util.o
	$(CC) -o $@ $^ $(LFLAGS)

fft.o: fft.c fft.h
	$(CC) -o $@ -c $(CFLAGS) $<

util.o: util.c util.h
	$(CC) -o $@ -c $(CFLAGS) $<

test_fft.o: test_fft.c fft.h util.h
	$(CC) -o $@ -c $(CFLAGS) $<

benchmark_fft.o: benchmark_fft.c fft.h util.h
	$(CC) -o $@ -c $(CFLAGS) $<

TestFFT: TestFFT.cpp fft.hpp util.hpp
	$(CXX) -o $@ $(CXXFLAGS) $<

TestFFTVector: TestFFTVector.cpp FFTVector.hpp fft.hpp util.hpp
	$(CXX) -o $@ $(CXXFLAGS) $<

BenchmarkFFT: BenchmarkFFT.cpp fft.hpp util.hpp
	$(CXX) -o $@ $(CXXFLAGS) $<

clean:
	rm -f *.o $(COBJS) $(CXXOBJS)
