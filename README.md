## Fast Fourier Transform

### C fast fourier transform library
### C++ fast fourier transform template library

#### `dft` Discrete Fourier Transform
- iterative
- N^2 time complexity

#### `idft` Inverse Discrete Fourier Transform
-  iterative
- N^2 time complexity

#### `fft_rec` Fast Fourier Transform
- recursive
- N log N time complexity
- based on the [Cooley-Tukey FFT algorithm] (Radix-2 DIT case)
    X[0,...,N−1] ← ditfft2(x, N, s):             DFT of (x[0], x[s], x[2s], ..., x[(N-1)s]):
        if N = 1 then
            X[0] ← x[0]                                    trivial size-1 DFT base case
        else
            X[0,...,N/2−1] ← ditfft2(x, N/2, 2s)           DFT of (x[0], x[2s], x[4s], ...)
            X[N/2,...,N−1] ← ditfft2(x+s, N/2, 2s)         DFT of (x[s], x[s+2s], x[s+4s], ...)
            for k = 0 to N/2−1 do                          combine DFTs of two halves into full DFT:
                t ← X[k]
                X[k] ← t + exp(−2πi k/N) X[k+N/2]
                X[k+N/2] ← t − exp(−2πi k/N) X[k+N/2]
            end for
        end if

#### `ifft_rec` Fast Fourier Transform
- recursive
- N log N time complexity

#### `fft_iter` Fast Fourier Transform
- iterative
- in place
- N log N time complexity
- based on [Cooley-Tukey FFT algorithm] (Data reordering, bit reversal, and in-place algorithms)
    algorithm iterative-fft is
        input: Array a of n complex values where n is a power of 2.
        output: Array A the DFT of a.

        bit-reverse-copy(a, A)
        n ← a.length
        for s = 1 to log(n) do
            m ← 2s
            ωm ← exp(−2πi/m)
            for k = 0 to n-1 by m do
                ω ← 1
                for j = 0 to m/2 – 1 do
                    t ← ω A[k + j + m/2]
                    u ← A[k + j]
                    A[k + j] ← u + t
                    A[k + j + m/2] ← u – t
                    ω ← ω ωm

        return A

#### `ifft_iter` Inverse Fast Fourier Transformn place
* iterative
* in place
* N log N time complexity

[Cooley-Tukey FFT algorithm]:<https://en.wikipedia.org/wiki/Cooley-Tukey_FFT_algorithm>
