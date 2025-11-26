#pragma once

#include <cmath>
#include <vector>
#include <array>
#include <functional>
#include <iomanip>
#include <iostream>

#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_linalg.h>

class Node {
public:
    unsigned int idx;
    double x;
    double y;
};

constexpr int nx = 3;
constexpr double x_min = 0;
constexpr double x_max = M_PI;
constexpr double dx = x_max / (nx - 1);

constexpr int ny = 3;
constexpr double y_min = 0;
constexpr double y_max = M_PI;
constexpr double dy = y_max / (ny - 1);

constexpr int M = (ny - 1) * (nx - 1); // liczba elementów
constexpr int num_nodes = nx*ny;       // liczba węzłów
constexpr int N_max = 4*nx*ny;         // liczba funkcji bazowych
constexpr double delta_x = 0.001;      // krok do obliczania pochodnych

constexpr int l_alpha_beta[4][2] = {{0, 0},{1, 0},{1, 1},{0, 1}};

extern std::array<Node,num_nodes> global_nodes;               // tablica węzłów globalnych
extern std::array<std::array<const Node*, 4>, M> local_nodes; // tablica węzłów lokalnych
extern std::vector<std::vector<double>> S;                    // globalna macierz sztywności
extern std::vector<double> F;                                 // globalny wektor obciążenia
extern std::vector<double> c;                                 // wektor rozwiazań 

// std::vector<std::vector<double>> A(n, std::vector<double>(n, 0.0));
// std::vector<double> c(n, 0.0); // wspolczynniki rozwiazania
// std::vector<double> b(n, 0.0);

/**
 * @brief Funkcja rho(x,y)
 */
double rho(double x, double y);

/**
 * @brief Dokladne rozwiazanie u(x,y)
 */
double u_analytical(double x, double y);

/**
 * @brief Oblicza iloraz roznicowy centralny pierwszego rzedu dla funkcji u w punkcie x
 */
double diff_u1(const std::function<double(double)>& u, double x, double dx = delta_x);

/**
 * @brief Oblicza iloraz roznicowy drugiego rzedu dla funkcji u w punkcie x
 */
double diff_u2(const std::function<double(double)>& u, double x, double dx = delta_x);

/**
 * @brief Rozwiazuje uklad rownan Ac = b metoda LU z biblioteki GSL
 */
void solve_Ac_b(std::vector<std::vector<double>>& A, std::vector<double>& c, std::vector<double>& b, int N);

/**
 * @brief Całkowanie numeryczne metodą kwadratury Gaussa-Legendre'a
 */
double integrate(std::function<double(double)> f, double a, double b, unsigned int n);

/**
 * @brief Funkcje ksztaltu Hermite'a dla elementu skonczonego trojpunktowego
 */
double shape_hermit_functions(int alpha, int i, double xi);

/**
 * @brief Funkcje wagowe dla elementu skonczonego czteropunktowego
 */
double weight_functions(int i, double xi_1, double xi_2);

/**
 * @brief Zwraca globalny indeks na podstawie indeksu lokalnego i indeksu funkcji kształtu
 */
int global_index(int local_index, int phi_i, int phi_j);

/**
 * @brief Tworzy tablice wezlow globalnych i lokalnych
 */
void make_nodes();

/**
 * @brief Polozenie globalne x zalezne od numeru elementu m i lokalnych współrzędnych xi_1, xi_2
 */
double x(int m, double xi_1, double xi_2);

/**
 * @brief Polozenie globalne y zalezne od numeru elementu m i lokalnych współrzędnych xi_1, xi_2
 */
double y(int m, double xi_1, double xi_2);

/**
 * @brief Składa globalna macierz sztywności i wektor obciążenia
 */
void assemble_global_matrices_MES_2D();

