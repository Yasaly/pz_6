#include <iostream>
#include <iomanip>
#include <vector>
#include <complex>
#include <chrono>
#include <fstream>
#include <cmath>
#include <windows.h>
#include "../headers/fourier.h"
#include <string>

#define PI 3.14159265358979323846

using namespace std;
using std::vector;
using std::complex;

string flag = "6";

void printSpectrumTable(const vector<complex<double>>& X,
                        const vector<complex<double>>& z)
{
    int N = static_cast<int>(X.size());
    int SETW = 22;

    cout << left << setw(SETW) << "m"
              << left << setw(SETW) << "Re z[m]"
              << left << setw(SETW) << "Re Z^[m]"
              << left << setw(SETW) << "Im Z^[m]"
              << left << setw(SETW) << "|Z^[m]|"
              << left << setw(SETW) << "phi"
              << '\n';

    cout << fixed << setprecision(6);

    for (int m = 0; m < N; ++m)
    {
        double amp   = abs(X[m]);
        double phase = arg(X[m]);

        cout << left << setw(SETW) << m
                  << left << setw(SETW) << z[m].real()
                  << left << setw(SETW) << X[m].real()
                  << left << setw(SETW) << X[m].imag()
                  << left << setw(SETW) << amp
                  << left << setw(SETW) << phase
                  << '\n';
    }
}

int main()
{
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);
    const int    N   = 512;
    const double A   = 3.0;
    const double B   = 0.26;
    const int    w1  = 1;
    const int    w2  = 193;
    const double phi = PI / 2.0;


    vector<complex<double>> z(N);

    //пункты 1-5
    if (flag == "1-5") {
        for (int j = 0; j < N; ++j)
        {
            double val = A * cos(2.0 * PI * w1 * j / N + phi)
                       + B * cos(2.0 * PI * w2 * j / N);
            z[j] = complex<double>(val, 0.0);
        }
    }

    //пункт 6
    else if (flag == "6") {
        for (int j = 0; j < N; ++j)
        {
            double val = 0.0;

            if (j >= N/4 && j < N/2)
            {
                val = A + B * cos(2.0 * PI * w2 * j / N);
            }
            else if (j >= 3*N/4 && j < N)
            {
                val = A + B * cos(2.0 * PI * w2 * j / N);
            }
            z[j] = complex<double>(val, 0.0);
        }
    }


    vector<complex<double>> Z_dft, Z_fft;

    auto t1 = chrono::high_resolution_clock::now();
    DFT(z, Z_dft);
    auto t2 = chrono::high_resolution_clock::now();

    auto t3 = chrono::high_resolution_clock::now();
    FFT(z, Z_fft);
    auto t4 = chrono::high_resolution_clock::now();

    auto dft_time = chrono::duration<double, milli>(t2 - t1).count();
    auto fft_time = chrono::duration<double, milli>(t4 - t3).count();

    cout << "Время:\n";
    cout << "DFT: " << dft_time << " мс\n";
    cout << "FFT: " << fft_time << " мс\n\n";

    printSpectrumTable(Z_dft, z);

    //Для получения чистого сигнала п.6
    // if (flag == "6") {
    //     int m_noise = w2;
    //     int band = 2;
    //
    //     for (int k = -band; k <= band; ++k)
    //     {
    //         int i1 = (m_noise + k + N) % N;
    //         int i2 = (N - m_noise + k + N) % N;
    //
    //         Z_fft[i1] = complex<double>(0.0, 0.0);
    //         Z_fft[i2] = complex<double>(0.0, 0.0);
    //     }
    //     IFFT(Z_fft, z);
    // }

    //Для получения чистого сигнала
    // Z_fft[193] = complex<double>(0.0, 0.0);
    // Z_fft[319] = complex<double>(0.0, 0.0);
    // IFFT(Z_fft, z);

    ofstream fout("signal.txt");
    if (fout)
    {
        fout << "# j  y\n";
        fout << setprecision(10);
        for (int j = 0; j < N; ++j)
            fout << j << " " << z[j].real() << "\n";

        fout.close();
        cout << "\nДанные записаны в 'signal.txt'.\n";
    }
    else
    {
        cerr << "Ошибка открытия файла 'signal.txt'.\n";
    }

    return 0;
}
