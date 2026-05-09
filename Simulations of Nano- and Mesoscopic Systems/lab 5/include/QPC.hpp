#pragma once
#include "params.hpp"
#include <string>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

/**
 * @brief Dyskretna siatka Nx x Ny z krokiem dx = dy
 */
struct Grid {
    double x_min, x_max;
    double y_min, y_max;
    int Nx, Ny;
    double m_eff;

    double dx() const { return (x_max - x_min) / (Nx + 1); }
    double dy() const { return (y_max - y_min) / (Ny + 1); }
    double alpha() const { return 0.5 / m_eff / (dy() * dy()); }
};

/**
 * @brief Mody propagujące się w prawo dla danej energii
 */
struct Modes {
    std::vector<Eigen::dcomplex> lambda; ///< wartości własne e^{ikx dx}
    std::vector<Eigen::VectorXcd> u;     ///< znormalizowane mody poprzeczne
    std::vector<double> v;               ///< prędkości v > 0
};

/**
 * @brief Rozwiązuje transport elektronowy w QPC metodą QTBM
 */
class QPCSolver {
public:
    explicit QPCSolver(Grid grid);

    /**
     * @brief Zad 1: relacja dyspersji E(kx) dla jednorodnego kanału
     */
    void dispersion(const std::string& filename, int N_kx = 500) const;

    /**
     * @brief Zad 2: mody propagujące się dla energii E [au]
     */
    Modes modes(double E, const std::string& filename   = "",
                          const std::string& filename_u = "") const;

    /**
     * @brief Zad 3: transmisja T i odbicie R metodą QTBM
     * @return {T, R}
     */
    std::pair<double, double> transmission(double E, double V_gates,
                                           const std::string& filename = "") const;

    /**
     * @brief Zad 4: mapa potencjału V(x,y) i T(Vgates)
     */
    void conductance_vs_gate(double E, double Vg_min, double Vg_max, int steps,
                             const std::string& filename,
                             const std::string& filename_V) const;

private:
    Grid m_grid;

    Eigen::MatrixXd make_tau() const;
    Eigen::MatrixXd make_H0(double E = 0.0) const;

    /**
     * @brief Potencjał QPC: V(x,y)
     */
    double V_gate(double x, double y, double V_gates) const;

    /**
     * @brief Macierz QTBM dla zadanej energii E i potencjału V_gates
     */
    Eigen::MatrixXcd make_QTBM(const Modes& m, double E, double V_gates = 0.0) const;

    std::pair<double, double> solve_TR(const Modes& m, const Eigen::MatrixXcd& M) const;
};
