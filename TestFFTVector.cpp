#include <iostream>
#include <vector>
#include "FFTVector.hpp"
#include "Util.hpp"

typedef std::vector<std::complex<double>> vec_type;

vec_type x = {1,2,3,4,3,2,1,0};

int main()
{
    util::print(std::cout, x, "x", 'f');

    vec_type X = fft::fft_iter(x);
    printf("X = fft_iter(x)...\n");
    util::print(std::cout, X, "X", 'f');

    vec_type x_ = fft::ifft_iter(X);
    printf("x = ifft_iter(X)...\n");
    util::print(std::cout, x_, "x", 'f');
}
