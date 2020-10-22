CC = gcc
CFLAGS = -std=c99 -pedantic-errors -Wall -O3
LFLAGS = -lm
OBJS = benchmark_fft.o test_fft.o fft.o util.o

all: test_fft benchmark_fft

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

clean:
	rm -f $(OBJS) benchmark_fft test_fft
