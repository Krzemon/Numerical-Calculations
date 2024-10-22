#include <iostream>
#include <fstream>
#include <cmath>
#include <random>    // dla generatora rozkładu normalnego
#include <string>

extern const double D;
extern const int N_max;
extern const double dt;
extern const double t_max;
extern const double sigma_dt;

extern double x[];
extern double y[];


inline double uniform() {
    return rand() / (double)RAND_MAX;
}

void gaus01(double& x, double& y);


// Funkcja symulująca krok dyfuzji
void SimStep();

// Funkcja zapisująca współczynniki do pliku
void saveCoefficients(double D_xx, double D_xy, double D_yy, int iteration);

// Funkcja zapisująca dane położenia do pliku
void saveData(const std::string& filename, const double x1[], const double y1[]);

// Funkcja obliczająca odległość między dwoma punktami
inline double distance(double x1, double y1, double x2, double y2) {
    return std::sqrt(std::pow(x2 - x1, 2) + std::pow(y2 - y1, 2));
}

// Funkcja obliczająca średnią arytmetyczną
double mean(const double x[], int n);

// Funkcja obliczająca średnią arytmetyczną kwadratów
double mean_square(const double x[], int n);

// Procedura symulująca proces dyfuzji z wykorzystaniem symulacji Wienera
void wiener();