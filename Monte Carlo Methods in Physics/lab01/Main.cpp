#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <ctime>

double uniform() {
    return (double)rand() / RAND_MAX;
}

double Bernoulli(double p) {
    if (uniform() <= p)
        return 1;
    else
        return 0;
}

void toFile(std::ofstream &file, int k, double relative_EX, double relative_std_dev) {
    const int DISPLAY_PRECISION = 7;
    file << std::fixed << std::setprecision(DISPLAY_PRECISION);
    file << k << '\t' << relative_EX << '\t' << relative_std_dev << std::endl;
}

void calculateAndSave(double p, std::ofstream &file) {
    const int N = 1e7;
    int k = 2;
    double x;
    double EX = 0.0;
    double EXsquare = 0.0;
    double relative_EX = 0.0;
    double relative_var = 0.0;

    srand(time(nullptr));

    for (int i = 1; i < N + 1; ++i) {
        x = Bernoulli(p);
        EX += x;
        EXsquare += (x * x);
        if (i == std::pow(10, k)) {
            double EX_temp = EX / i;
            double EXsquare_temp = EXsquare / i;
            double variance = EXsquare_temp - EX_temp * EX_temp;

            relative_EX = std::fabs((p - EX_temp) / p);
            relative_var = std::fabs((variance - (p*(1-p))) / (p*(1-p)));
            toFile(file, k, relative_EX, relative_var);
            ++k;
        }
    }
}

int main() {
    double tab[] = {0.1, 0.5, 0.9};

    std::ofstream file_p01("../data/p01.dat", std::ios::out);
    std::ofstream file_p05("../data/p05.dat", std::ios::out);
    std::ofstream file_p09("../data/p09.dat", std::ios::out);

    calculateAndSave(tab[0], file_p01);
    calculateAndSave(tab[1], file_p05);
    calculateAndSave(tab[2], file_p09);

    file_p01.close();
    file_p05.close();
    file_p09.close();

    return 0;
}
