#include "../headers/fourier.h"
#include <cmath>

#define PI 3.14159265358979323846

using std::complex;
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


void FFT(const vector<complex<double>>& x,
         vector<complex<double>>& X)
{
    int N = static_cast<int>(x.size());
    X.assign(N, complex<double>(0.0, 0.0));

    if (N == 1)
    {
        X[0] = x[0];
        return;
    }

    int M = N / 2;
    vector<complex<double>> x_u(M), x_v(M);
    for (int i = 0; i < M; ++i)
    {
        x_u[i] = x[2 * i];
        x_v[i]  = x[2 * i + 1];
    }

    vector<complex<double>> X_u, X_v;
    FFT(x_u, X_u);
    FFT(x_v,  X_v);

    for (int m = 0; m < M; ++m)
    {
        double angle = -2.0 * PI * m / N;
        complex<double> w(cos(angle), sin(angle));
        complex<double> t = w * X_v[m];

        X[m]     = X_u[m] + t;
        X[m + M] = X_u[m] - t;
    }
}

void IFFT(const vector<complex<double>>& X,
          vector<complex<double>>& x)
{
    int N = static_cast<int>(X.size());
    vector<complex<double>> X_conj(N);

    for (int i = 0; i < N; ++i)
        X_conj[i] = conj(X[i]);

    vector<complex<double>> tmp;
    FFT(X_conj, tmp);

    x.assign(N, complex<double>(0.0, 0.0));
    for (int i = 0; i < N; ++i)
        x[i] = conj(tmp[i]) / static_cast<double>(N);
}
