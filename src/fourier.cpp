#include "../headers/fourier.h"
#include <cmath>
#include <complex>

#define PI 3.14159265358979323846

using namespace std;

void DFT(const vector<complex<double>>& x,
         vector<complex<double>>& X)
{
    int N = static_cast<int>(x.size());
    X.assign(N, complex<double>(0.0, 0.0));

    for (int k = 0; k < N; ++k)
    {
        complex<double> sum(0.0, 0.0);
        for (int n = 0; n < N; ++n)
        {
            double angle = -2.0 * PI * k * n / N;
            complex<double> w(cos(angle), sin(angle));
            sum += x[n] * w;
        }
        X[k] = sum;
    }
}

void IDFT(const vector<complex<double>>& X,
          vector<complex<double>>& x)
{
    int N = static_cast<int>(X.size());
    x.assign(N, complex<double>(0.0, 0.0));

    for (int n = 0; n < N; ++n)
    {
        complex<double> sum(0.0, 0.0);
        for (int k = 0; k < N; ++k)
        {
            double angle = 2.0 * PI * k * n / N;
            complex<double> w(cos(angle), sin(angle));
            sum += X[k] * w;
        }
        x[n] = sum / static_cast<double>(N);
    }
}


void FFT(vector<complex<double>>& x, int N, int start, int step){
    if (N == 1) return;
    int M = N / 2;

    FFT(x, M, start, step * 2);
    FFT(x, M, start + step, step * 2);
    double angle_step = -2.0 * PI / static_cast<double>(N);
    complex<double> w_step(cos(angle_step), sin(angle_step));
    complex<double> w(1.0, 0.0);

    for (int m = 0; m < M; ++m)
    {
        int idx_u = start + step * (2 * m);
        int idx_v = start + step * (2 * m + 1);
        complex<double> u = x[idx_u];
        complex<double> v = x[idx_v];
        complex<double> t = w * v;
        x[idx_u] = u + t;
        x[idx_v] = u - t;
        w *= w_step;
    }
}

void IFFT(const vector<complex<double>>& X,
          vector<complex<double>>& x)
{
    x = X;
    int N = static_cast<int>(x.size());
    if (N == 0) return;

    for (int i = 0; i < N; ++i)
        x[i] = conj(x[i]);

    FFT(x, N, 0, 1);

    for (int i = 0; i < N; ++i)
        x[i] = conj(x[i]) / static_cast<double>(N);
}
