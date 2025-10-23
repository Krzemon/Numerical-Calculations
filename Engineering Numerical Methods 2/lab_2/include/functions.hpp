#pragma once

#include <iostream>
#include <iomanip>
#include <cmath>
#include <functional>
#include <vector>
#include <array>

#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_linalg.h>

constexpr double delta_x = 0.001;

constexpr auto f = [](double x) -> double {
    return 2*(1 - 3*x + x*x - M_PI*M_PI) * std::exp(-x);
};

/**
 * @brief Oblicza rozwiazanie analityczne w punkcie x
 */
inline double analytical_solution(double x) {
    return (x - M_PI)*(x + M_PI)*std::exp(-x);
}

/**
 * @brief Oblicza iloraz roznicowy centralny pierwszego rzedu dla funkcji u w punkcie x
 */
inline double diff_u1(std::function<double(double)>& u, double x, double dx = delta_x) {
    return (u(x + delta_x) - u(x - delta_x)) / (2*delta_x);
}

/**
 * @brief Oblicza iloraz roznicowy drugiego rzedu dla funkcji u w punkcie x
 */
inline double diff_u2(std::function<double(double)>& u, double x, double dx = delta_x) {
    return (u(x + delta_x) - 2*u(x) + u(x - delta_x)) / (delta_x*delta_x);
}

/**
 * @brief Oblicza funkcje bazowa 1 dla indeksu i w punkcie x
 */
inline auto gen_base_1(int N) -> std::vector<std::function<double(double)>> {
    std::vector<std::function<double(double)>> v;
    for (int i = 1; i <= N; ++i) {
        v.push_back([i](double x) {
            return std::cos((i - 0.5)*x) * std::exp(-x);
        });
    }
    return v;
}

/**
 * @brief Oblicza funkcje bazowa 2 dla indeksu i w punkcie x
 */
inline auto gen_base_2(int N) -> std::vector<std::function<double(double)>> {
    std::vector<std::function<double(double)>> v;
    for (int i = 1; i <= N; ++i) {
        v.push_back([i](double x) {
            return (x - M_PI)*(x + M_PI)*std::pow(x, i-1);
        });
    }
    return v;
}

/**
 * @brief Sprawdza warunki brzegowe dla funkcji bazowej w punkcie x = -pi oraz x = pi
 */
inline void check_bc(const std::vector<std::function<double(double)>>& v, int N, int base_type) {
    std::vector<double> v_vec_a(N), v_vec_b(N);
    double sum_a = 0.0, sum_b = 0.0;
    for (int i = 0; i < N; ++i) {
        v_vec_a[i] = v[i](-M_PI);
        v_vec_b[i] = v[i](M_PI);
        sum_a += v_vec_a[i];
        sum_b += v_vec_b[i];
    }
    std::cout << "Wartosc funkcji dla bazy " << base_type << " w x = -pi: " << sum_a << std::endl;
    std::cout << "Wartosc funkcji dla bazy " << base_type << " w x =  pi: " << sum_b << std::endl;
}

/**
 * @brief Sprawdza maksymalne residuum dla ukladu rownan Ac = b
 */
inline double check_max_row_residuum(const std::vector<std::vector<double>>& A, const std::vector<double>& c, const std::vector<double>& b, int n) {
    double max_abs = 0;
    for (int i = 0; i < n; ++i){
        double s = 0;
        for (int j = 0; j < n; ++j) 
            s += A[i][j] * c[j];
        double ri = s - b[i];
        max_abs = std::max(max_abs, std::abs(ri));
    }
    return max_abs;
}

/**
 * @brief Rozwiazuje uklad rownan Ac = b metoda LU z biblioteki GSL
 */
inline void solve_Ax_b(std::vector<std::vector<double>>& A, std::vector<double>& c, std::vector<double>& b, int N) {
    gsl_matrix *m = gsl_matrix_alloc(N, N);
    gsl_vector *B = gsl_vector_alloc(N);
    gsl_vector *x = gsl_vector_alloc(N);

    for (int i = 0; i < N; ++i) {
        gsl_vector_set(B, i, b[i]);
        for (int j = 0; j < N; ++j)
            gsl_matrix_set(m, i, j, A[i][j]);
    }
    // Rozwiązywanie (LU)
    gsl_permutation *p = gsl_permutation_alloc(N);
    int signum;
    gsl_linalg_LU_decomp(m, p, &signum);
    gsl_linalg_LU_solve(m, p, B, x);

    for (int i = 0; i < N; ++i)
        c[i] = gsl_vector_get(x, i);

    gsl_matrix_free(m);
    gsl_vector_free(B);
    gsl_vector_free(x);
    gsl_permutation_free(p);
}

/**
 * @brief Generuje wektor x od l_edge do r_edge z krokiem dx
 */
inline auto gen_x(double l_edge, double r_edge, double dx = delta_x) -> std::vector<double> {
    std::vector<double> x_vec;
    while (l_edge <= r_edge) {
        x_vec.push_back(l_edge);
        l_edge += dx;
    }
    return x_vec;
}