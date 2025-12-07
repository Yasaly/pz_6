#ifndef FOURIER_H
#define FOURIER_H

#include <vector>
#include <complex>

void DFT(const std::vector<std::complex<double>>& x,
         std::vector<std::complex<double>>& X);

void IDFT(const std::vector<std::complex<double>>& X,
          std::vector<std::complex<double>>& x);

void FFT(const std::vector<std::complex<double>>& x,
         std::vector<std::complex<double>>& X);

void IFFT(const std::vector<std::complex<double>>& X,
          std::vector<std::complex<double>>& x);

#endif // FOURIER_H
