#pragma once

#include <cmath>
#include <vector>
#include <iomanip>
#include <iostream>

#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_linalg.h>

extern double L;

extern double x_min;
extern double x_max;
extern int nx;
extern double dx;

extern double y_min;
extern double y_max;
extern int ny;
extern double dy;

struct Node {
    int index;   // globalny indeks
    double x,y;  // polozenie wezla
};
struct Element { 
    int i, j, k; // globalne indeksy wezlow elementu
    int a, b, c; // lokalne  indeksy wezlow elementu
    int index;   // indeks elementu
};

extern std::vector<Node> nodes;
extern std::vector<Element> elements;

/**
 * @brief Funkcja rho(x,y) - rozklad ladunku 
 */
double rho(double x, double y);

/**
 * @brief Rozwiazuje uklad rownan Ac = b metoda LU z biblioteki GSL
 */
void solve_Ac_b(std::vector<std::vector<double>>& A, std::vector<double>& c, std::vector<double>& b, int N);

/**
 * @brief Tworzy siatke zlozona z wezlow i elementow
 */
void make_grid(int new_nx = nx, int new_ny = ny, double new_x_min = x_min, double new_x_max = x_max, 
                                                 double new_y_min = y_min, double new_y_max = y_max);

/**
 * @brief Zwraca wartosc liniowej funkcji ksztaltu w przestrzeni referencyjnej
 */
void shape_functions(double z, double e, double &phi0, double &phi1, double &phi2);

/**
 * @brief Oblicza jakobian przeksztalcenia
 */
void jacobian(const Node &A, const Node &B, const Node &C, double dzeta, double eta, double (&J)[2][2], double (&invJ)[2][2], double &detJ);

/**
 * @brief Sklada lokalna macierz sztywnosci dla elementu trojkatnego
 */
void assemble_local_matrice_triangle(const std::vector<Node> &nodes, const Element &el, double (&local_matrix)[3][3]);

/**
 * @brief Sklada lokalny wektor obciazen dla elementu trojkatnego
 */
void assemble_local_vector_triangle(const std::vector<Node> &nodes, const Element &el, double (&local_vector)[3]);

/**
 * @brief Zwraca globalny numer wezla na podstawie lokalnego numeru wezla oraz numeru elementu
 */
int nr(int local_index, int element_index);

// bool isBoundaryNode(const Node &p);