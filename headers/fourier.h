#ifndef FOURIER_H
#define FOURIER_H

#include <vector>
#include <complex>

using namespace std;

void DFT(const vector<complex<double>>& x,
         vector<complex<double>>& X);

void IDFT(const vector<complex<double>>& X,
          vector<complex<double>>& x);

void FFT(vector<complex<double>>& x, int N, int start, int step);

void IFFT(const vector<complex<double>>& X,
          vector<complex<double>>& x);

#endif // FOURIER_H
