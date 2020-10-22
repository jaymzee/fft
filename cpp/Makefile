CC = gcc
CXX = g++
CFLAGS = -std=c++11 -pedantic-errors -Wall -O3
LFLAGS =
OBJS = TestFFT TestFFTVector BenchmarkFFT

all: $(OBJS)

TestFFT: TestFFT.cpp fft.hpp util.hpp
	$(CXX) -o $@ $(CFLAGS) $<

TestFFTVector: TestFFTVector.cpp FFTVector.hpp fft.hpp util.hpp
	$(CXX) -o $@ $(CFLAGS) $<

BenchmarkFFT: BenchmarkFFT.cpp fft.hpp util.hpp
	$(CXX) -o $@ $(CFLAGS) $<

clean:
	rm -f $(OBJS) *.o
