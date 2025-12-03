#pragma once
#include <vector>
#include <array>
#include <functional>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <iomanip>

#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>

extern "C" {
    void dsygvd_(int* itype, char* jobz, char* uplo, int* n,
                 double* a, int* lda, double* b, int* ldb,
                 double* w, double* work, int* lwork, int* iwork, int* liwork, int* info);
}

struct Node {
    double x,y;
    int idx;
};

using NodeSet = std::array<const Node*, 3>;

extern int g_nx, g_ny;
extern double g_xmin, g_xmax, g_ymin, g_ymax;
extern double g_dx, g_dy;
extern int g_N; // liczba węzłów
extern int g_M; // liczba elementów

extern std::vector<Node> g_nodes;
extern std::vector<NodeSet> g_elem_nodes;
extern std::vector<double> g_E; // macierz globalna sztywności
extern std::vector<double> g_O; // macierz globalna przekrywania
extern std::vector<std::array<double,9>> g_E_local; // macierz lokalna sztywności
extern std::vector<std::array<double,9>> g_O_local; // macierz lokalna przekrywania

extern std::vector<double> evals;                   // wartosci wlasne
extern std::vector<std::vector<double>> evecs;      // wektory wlasne

constexpr double DXI = 0.001;
/**
 * @brief Funkcje ksztaltu w przestrzeni referencyjnej
 */
static const std::function<double(double,double)> phi[3] = {
    [](double dzeta, double eta){ return -0.5*(dzeta + eta); },
    [](double dzeta, double eta){ return  0.5*(1 + dzeta);   },
    [](double dzeta, double eta){ return  0.5*(1 + eta);     }
};

/**
 * @brief Oblicza pochodna funkcji f(x) w punkcie x
 */
inline double df_num(std::function<double(double)> f, double x) {
    double h = DXI;
    return (f(x + h) - f(x - h)) / (2.0*h);
}

/**
 * @brief Oblicza pochodna dwuwymiarowa funkcji f(x,y) w kierunku x w punkcie x,y
 */
inline double df_num_2D_x(std::function<double(double,double)> f, double x, double y) {
    double h = DXI;
    return (f(x + h, y) - f(x - h, y)) / (2.0*h);
}

/**
 * @brief Oblicza pochodna dwuwymiarowa funkcji f(x,y) w kierunku y w punkcie x,y
 */
inline double df_num_2D_y(std::function<double(double,double)> f, double x, double y) {
    double h = DXI;
    return (f(x, y + h) - f(x, y - h)) / (2.0*h);
}
/**
 * @brief Zwraca globalny indeks wezla na podstawie indeksu elementu oraz lokalnego indeksu wezla
 */
const int lg(int elem_idx, int local_node_idx); 

/**
 * @brief Rozwiazuje uklad rownan Ac = b metoda LU z biblioteki GSL
 */
// void solve_Ac_b(std::vector<std::vector<double>>& A, std::vector<double>& c, std::vector<double>& b, int N);

/**
 * @brief Calkowanie w przestrzeni referencyjnej
 */
double integrate_over_reference_space(std::function<double(double,double)> f);

/**
 * @brief Mapuje element na punkt w przestrzeni fizycznej x oraz y
 */
double x_map(int elem_idx, double xi1, double xi2);
double y_map(int elem_idx, double xi1, double xi2);

/**
 * @brief Oblicza jakobian przeksztalcenia
 */
double jacobian(int elem_idx, double xi1, double xi2);

/**
 * @brief Tworzy siatke zlozona z wezlow i elementow
 */
void make_grid(int nx, double x_min, double x_max,
               int ny, double y_min, double y_max);

/**
 * @brief Wypelnia dane wezlow
 */
void fill_node_data();

/**
 * @brief Wypelnia lokalne macierze
 */
void fill_local_matrices();

/**
 * @brief Sklada lokalna macierz sztywnosci dla elementu trojkatnego
 */
void assemble_matrices();

/**
 * @brief Uwzglednia warunki brzegowe
 */
void apply_boundary_conditions();

/**
 * @brief Calkowanie metoda newmarka
 */
// void newmark_time_integration(double dt, int n_steps, const std::vector<double>& u0, std::vector<double>& u_curr);

/**
 * @brief Rozwiazuje problem wlasny z wykorzystaniem lapack
 */
void solve_generalized_eigen_lapack(const double* E, const double* O, int N, std::vector<double>& evals, std::vector<std::vector<double>>& evecs);



std::vector<double> solve_linear_system(const std::vector<double>& A, const std::vector<double>& b);