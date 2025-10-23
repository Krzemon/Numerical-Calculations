#include "functions.hpp"
#include <iostream>
#include <iomanip>
#include <array>
#include <vector>
#include <cmath>

constexpr double dE = 0.2;
constexpr double E_max = 150.0;

// ------------------ Zadanie 1 ------------------
inline void calc_zad_1(std::ostream& file, const int& l = 0) {
    std::array<double, N+1> R{};
    for (double E = dE; E <= E_max+dE/10; E += dE) {
        calc_Rnum(l, delta_r, E, R.data());
        file << std::fixed << std::setprecision(6) << E << " " << R[N] << "\n";
    }
}

// ------------------ Zadanie 2 ------------------
inline void calc_zad_2(std::ostream& file, const int& l = 0) {
    std::array<double, N+1> R{};
    file << "E_num\t\tE_theor\n";

    std::cout << "------------------ZADANIE 2------------------\n";
    std::cout << "------------------  l = "<< l << " ------------------\n";
    std::cout << "E_num\t\tE_theor\n";
    
    double E_prev = dE;
    calc_Rnum(l, delta_r, E_prev, R.data());
    double R_prev = R[N];

    std::vector<double> E_num_found;

    for (double E = 2*dE; E <= E_max && E_num_found.size()<p_size; E+=dE) {
        calc_Rnum(l, delta_r, E, R.data());
        double R_curr = R[N];

        if (R_prev * R_curr < 0.0) {
            // metoda siecznych
            double E1 = E_prev, E2 = E;
            double R1 = R_prev, R2 = R_curr;
            double tol = 1e-6;
            int iter = 0, max_iter = 100;
            double E_next = E2;
            while (std::abs(E2-E1) > tol && iter<max_iter) {
                E_next = E2 - R2*(E2-E1)/(R2-R1);
                calc_Rnum(l, delta_r, E_next, R.data());
                double R_next = R[N];
                E1 = E2; R1 = R2;
                E2 = E_next; R2 = R_next;
                iter++;
            }
            E_num_found.push_back(E_next);
            double E_theor = 0.5*std::pow(alpha_0[E_num_found.size()-1]/L,2.0);\
            file << std::fixed << std::setprecision(6) << E_next << "\t"<< E_theor << "\n";
            std::cout << std::fixed << std::setprecision(6) << E_next << "\t"<< E_theor << "\n";
        }
        E_prev = E;
        R_prev = R_curr;
    }
    std::cout << "---------------------------------------------\n";
}

// ------------------ Zadanie 3 ------------------
inline void calc_zad_3(std::ostream& file, const int& l = 0) {
    std::array<double, p_size> E_theor{};
    for (size_t i = 0; i < p_size; ++i)
        E_theor[i] = calc_E(l, i+1, L);

    std::array<std::array<double, n_size()>, p_size> deltaR{};

    for (size_t p = 0; p < p_size; ++p) {
        std::array<double, n_size()> R_num{};
        std::array<double, n_size()> R_exact{};

        calc_Rnum(l, delta_r, E_theor[p], R_num.data());
        for (size_t i = 0; i < n_size(); ++i)
            R_exact[i] = calc_R(l, p+1, r_nodes[i], L);

        for (size_t i = 0; i < n_size(); ++i)
            deltaR[p][i] = R_exact[i] - R_num[i];
    }

    file << "r";
    for (size_t p = 1; p <= p_size; ++p)
        file << "\tDeltaR_p" << p;
    file << "\n";

    for (size_t i = 0; i < n_size(); ++i) {
        file << std::fixed << std::setprecision(6) << r_nodes[i];
        for (size_t p = 0; p < p_size; ++p)
            file << "\t" << std::fixed << std::setprecision(6) << deltaR[p][i];
        file << "\n";
    }
}

// ------------------ Zadanie 4 ------------------

inline void calc_zad_4_1(std::ostream& file, const int& l = 1) {
    std::vector<double> U(n_size(), 0.0);
    std::vector<double> R(n_size(), 0.0);

    file << std::fixed << std::setprecision(6);
    for (double E = dE; E <= E_max+dE/10; E += dE) {
        calc_U_Numerov(l, E, U);
        normalize_R(U, R, l);
        file << E << " " << R.back() << "\n";
    }
}

inline void calc_zad_4_2(const int& l = 1) {
    std::vector<double> U(n_size(), 0.0);
    std::vector<double> R(n_size(), 0.0);
    
    std::cout << "------------------ZADANIE 4------------------\n";
    std::cout << "------------------  l = "<< l << " ------------------\n";
    std::cout << "E_num\t\tE_theor\n";
    
    double E_prev = dE;
    calc_U_Numerov(l, E_prev, U);
    normalize_R(U, R, l);
    double R_prev = R.back();

    std::vector<double> E_num_found;

    for (double E = 2 * dE; E <= E_max && E_num_found.size() < p_size; E += dE) {
        calc_U_Numerov(l, E, U);
        normalize_R(U, R, l);
        double R_curr = R.back();

        if (R_prev * R_curr < 0.0) {
            // metoda siecznych
            double E1 = E_prev, E2 = E;
            double R1 = R_prev, R2 = R_curr;
            double tol = 1e-6;
            int iter = 0, max_iter = 100;
            double E_next = E2;

            while (std::abs(E2 - E1) > tol && iter < max_iter) {
                E_next = E2 - R2 * (E2 - E1) / (R2 - R1);
                calc_U_Numerov(l, E_next, U);
                normalize_R(U, R, l);
                double R_next = R.back();
                E1 = E2; R1 = R2;
                E2 = E_next; R2 = R_next;
                iter++;
            }
            E_num_found.push_back(E_next);
        }

        E_prev = E;
        R_prev = R_curr;
    }

    for (size_t p = 0; p < E_num_found.size(); ++p) {
        double E_theor = calc_E(l, p + 1, L);
        std::cout << std::fixed << std::setprecision(6)
                  << E_num_found[p] << "\t" << E_theor << "\n";
    }
    std::cout << "---------------------------------------------\n";
}

inline void calc_zad_4_3(std::ostream& file, const int& l = 1) {
    std::vector<double> U(n_size(), 0.0);
    std::vector<double> R_num(n_size(), 0.0);
    std::vector<double> R_exact(n_size(), 0.0);

    const int n_plot = 4;
    std::vector<std::vector<double>> deltaR_all(n_plot, std::vector<double>(n_size(), 0.0));

    file << "r";
    for (int p = 1; p <= n_plot; ++p)
        file << "\tDeltaR_p" << p;
    file << "\n";

    for (int p = 0; p < n_plot; ++p) {
        double E = calc_E(l, p + 1, L);

        calc_U_Numerov(l, E, U);
        normalize_R(U, R_num, l);

        for (size_t i = 0; i < n_size(); ++i)
            R_exact[i] = calc_R(l, p + 1, r_nodes[i], L);

        for (size_t i = 0; i < n_size(); ++i)
            deltaR_all[p][i] = R_exact[i] - R_num[i];
    }

    for (size_t i = 0; i < n_size(); ++i) {
        file << std::fixed << std::setprecision(6) << r_nodes[i];
        for (int p = 0; p < n_plot; ++p)
            file << "\t" << std::fixed << std::setprecision(6) << deltaR_all[p][i];
        file << "\n";
    }
}
