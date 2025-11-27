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

constexpr double x_min = 0;
constexpr double x_max = M_PI;
extern int nx;
extern double dx;

constexpr double y_min = 0;
constexpr double y_max = M_PI;
extern int ny;
extern double dy;

extern int M;           // liczba elementów
extern int num_nodes;   // liczba węzłów
extern int N_max;       // liczba funkcji bazowych
extern double delta_x;  // krok do obliczania pochodnych

constexpr int l_alpha_beta[4][2] = {{0, 0},{1, 0},{1, 1},{0, 1}};

extern std::vector<Node> global_nodes;                        // tablica węzłów globalnych
extern std::vector<std::array<const Node*, 4>> local_nodes;   // tablica węzłów lokalnych

extern std::vector<std::vector<double>> S;                    // globalna macierz sztywności
extern std::vector<double> F;                                 // globalny wektor obciążenia
extern std::vector<double> c;                                 // wektor rozwiazań 

extern std::vector<std::vector<double>> S_original;
extern std::vector<double> F_original;

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
 * @brief Aktualizuje parametry siatki
 */
void update_params(int new_nx, int new_ny);

/**
 * @brief Składa globalna macierz sztywności i wektor obciążenia
 */
void assemble_global_matrices_MES_2D();

/**
 * @brief Modyfikuje macierz sztywnosci S oraz wektor obciazen F o warunki brzegowe
 */
void border_conditions();

/**
 * @brief Oblicza calke funkcjonalna (Zasada Rayleigha-Ritz'a)
 */
double functional_integral();

/**
 * @brief Rozwiązanie numeryczne u(x,y) na podstawie wektora współczynników c
 */
double u_solution(double x, double y);