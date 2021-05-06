#include <iostream>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <fstream>
#include <complex>
#include <valarray>
#include <iomanip>
#include <pybind11/pybind11.h>
#include<pybind11/complex.h>
#include <pybind11/stl.h>

double pi = 3.14159265359;
const double PI = 3.141592653589793238460;
using cd = std::complex<long double>;



std::vector<cd> FFT(std::vector<cd> &fftArray)
{
    int N = fftArray.size();
    if (N == 1)
        return fftArray;
    std::vector<cd> even(N/2,0);
    std::vector<cd> odd(N/2,0);
    for (int i = 0; i != N / 2; i++)
    {
        even[i] = fftArray[2 * i];
        // even.push_back(x);
        odd[i] = fftArray[2 * i + 1];
        // odd.push_back(x);
    }

    // conquer
    std::vector<cd>Feven(N/2,0);
    Feven = FFT(even);
    std::vector<cd>Fodd(N/2,0);
    Fodd = FFT(odd);

    // combine
    std::vector<cd>freqs(N,0);
    for (size_t k = 0; k < N / 2; ++k)
    {
        cd t = (cd(cos(-2 * PI * double(k) / N), sin(-2 * PI * double(k) / N))) * Fodd[k];
        freqs[k] = Feven[k] + t;
        freqs[k + N / 2] = Feven[k] - t;
    }
    return freqs;
}

std::vector<cd> DFT(std::vector<double> signal)
{
    std::vector<cd> ftArray;
    int N = signal.size();
    for (int k = 0; k < N; k++)
    {
        long double reals = 0.0;
        long double imags = 0.0;
        for (int j = 0; j < N; j++)
        {
            long double angleTerm = 2 * pi * j * k * (1.0 / N);
            long double cosine = cos(angleTerm);
            long double sine = sin(angleTerm);
            reals += (signal[j] * cosine);
            imags += (signal[j] * -1 * sine);
        }
        // std::cout<<imags<<std::endl;
        cd sum;
        sum = cd(reals, imags);
        ftArray.push_back(sum);
    }
    return ftArray;
}

// int main()
// {
//     std::vector<double> signal;
//     std::ifstream inFile;
//     inFile.open("wave.txt");
//     if (!inFile)
//     {
//         std::cout << "Unable to open the file";
//         exit(1);
//     }
//     double dat;
//     while (inFile >> dat)
//     {
//         double data = dat;
//         signal.push_back(data);
//     }
//     std::vector<cd> vfft;
//     std::vector<cd> vdft;
//     for (int i = 0; i < signal.size(); i++)
//     {
//         cd var = signal[i];
//         vfft.push_back(var);
//     }

//     vdft = DFT(signal);
//     FFT(vfft);

//     std::cout << "fft" << std::endl;
//     for (int i = 0; i < 8; ++i)
//     {
//         std::cout << std::setprecision(9) << signal[i] << std::endl;
//     }
//     std::cout << signal.size() << std::endl;
// }

PYBIND11_MODULE(fourier, m) {
    // m.doc() = "pybind11 example plugin"; // optional module docstring
    m.def("DFT", &DFT, "Fourier transform");
    m.def("FFT", &FFT, "Fast fourier transform");
}