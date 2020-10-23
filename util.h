#ifndef FFT_UTIL_H_INCLUDED
#define FFT_UTIL_H_INCLUDED

#include <complex.h>

void print_complex(double complex x, char format);

void print_complex_array(double complex *x, int length,
                         const char *name, char fmt);

#endif
