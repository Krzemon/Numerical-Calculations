#pragma once

#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <array>
#include <iomanip>
#include <functional>
#include <gsl/gsl_sf_bessel.h>

constexpr double L = 1.0;
constexpr int N = 100;
constexpr double delta_r = L / N;
constexpr double alpha_0[] = {2.4048, 5.5200, 8.6537, 11.7915};
constexpr double alpha_1[] = {3.8317, 7.0155, 10.1734, 13.3236};
constexpr size_t p_size = sizeof(alpha_0)/sizeof(alpha_0[0]);

constexpr size_t n_size() {
    return N+1;
}
constexpr auto r_nodes = [] {
    std::array<double, n_size()> arr{};
    for (size_t i = 0; i < n_size(); ++i)
        arr[i] = i * delta_r;
    return arr;
}();

/**
 * @brief Oblicza wartosc funkcji R_l(r) = J_|l|(α_{l,p} * r / L)
 * @param l moment pedu
 * @param p numer zera funkcji Bessela
 * @param r odleglosc
 * @param L dlugosc
 */
inline double calc_R(int l, int p, double r, double L) {
    double alpha = (l == 0) ? alpha_0[p - 1] : alpha_1[p - 1];
    double x = alpha * r / L;
    if (l == 0)
        return gsl_sf_bessel_J0(x);
    else if (l == 1)
        return gsl_sf_bessel_J1(x);
    else
        return gsl_sf_bessel_Jn(l, x);
}

/**
 * @brief Oblicza wartosc energii E_{l,p}
 * @param l moment pedu
 * @param p numer zera funkcji Bessela
 * @param L dlugosc
 */
inline double calc_E(int l, int p, double L) {
    double alpha = (l == 0) ? alpha_0[p - 1] : alpha_1[p - 1];
    return 0.5 * pow(alpha / L, 2);
}

/**
 * @brief Oblicza numeryczne R_n(r) dla zadanej energii E metodą roznic skonczonych
 * @param l moment pedu
 * @param delta_r krok siatki
 * @param E energia
 * @param R tablica wynikowa R[0..N]
 */
inline void calc_Rnum(int l, double delta_r, double E, double R[]) {
    if (delta_r == 0.0) {
        std::cerr << "Error: delta_r nie moze byc zerem.\n";
        return;
    }
    R[0] = 1.0;
    R[1] = 1.0;
    for (int i = 1; i < N; ++i) {
        double ri = i * delta_r;
        double A = 1.0/(delta_r*delta_r) + 1.0/(2.0*ri*delta_r);
        double B = 2.0/(delta_r*delta_r) + (l*l)/(ri*ri) - 2.0*E;
        double C = -1.0/(delta_r*delta_r) + 1.0/(2.0*ri*delta_r);
        R[i+1] = (B*R[i] + C*R[i-1])/A;
    }
}

/**
 * @brief Oblicza numeryczne U(r) metodą Numerova
 * @param l moment pedu
 * @param delta_r krok siatki
 * @param E energia
 * @param U wektor wynikowy U[0..N]
 * @param N liczba punktow siatki
 */
inline void calc_U_Numerov(int l, double E, std::vector<double>& U) {
    if (n_size() < 3) return;

    U[0] = 0.0;
    U[1] = 1.0;
    auto g = [&](double r) {
        return (r != 0.0) ? (1 - 4*l*l) / (4 * r * r) + 2.0 * E : 2.0 * E;
    };

    for (size_t i = 1; i < n_size()-1; ++i) { // ... N --> U[i+1] = U[N]
        double g_im1 = g(r_nodes[i - 1]);
        double g_i   = g(r_nodes[i]);
        double g_ip1 = g(r_nodes[i + 1]);

        U[i + 1] = (2.0 * (1.0 - 5.0 * delta_r * delta_r / 12.0 * g_i) * U[i]
                    - (1.0 + delta_r * delta_r / 12.0 * g_im1) * U[i - 1])
                   / (1.0 + delta_r * delta_r / 12.0 * g_ip1);

        if (std::isnan(U[i + 1]) || std::isinf(U[i + 1]))
            U[i + 1] = 0.0;
    }
}

/**
 * @brief Normalizuje funkcje falowa R(r) na podstawie U(r)
 * @param U wektor U[0..N]
 * @param R wektor wynikowy R[0..N]
 * @param delta_r krok siatki
 * @param N liczba punktow siatki
 * @param l moment pedu
 */
inline void normalize_R(const std::vector<double>& U, std::vector<double>& R, int l=1) {
    if (l == 0){
        R[0] = 1.0 / std::sqrt(delta_r);
        R[1] = 1.0 / std::sqrt(delta_r);
        for (size_t i = 2; i < n_size(); ++i) {
            double r = r_nodes[i];
            R[i] = (r != 0.0) ? U[i] / std::sqrt(r) : 0.0;
        }
    } else{
        R[0] = 0.0;
        for (size_t i = 2; i < n_size(); ++i) {
            double r = r_nodes[i];
            R[i] = (r != 0.0) ? U[i] / std::sqrt(r) : 0.0;
        }
    }

    double norm = 0.0;
    for (size_t i = 0; i < n_size(); ++i)
        norm += R[i] * R[i] * r_nodes[i] * delta_r;

    norm = std::sqrt(norm);
    if (norm == 0.0) norm = 1.0;

    for (size_t i = 0; i < n_size(); ++i)
        R[i] /= norm;
}