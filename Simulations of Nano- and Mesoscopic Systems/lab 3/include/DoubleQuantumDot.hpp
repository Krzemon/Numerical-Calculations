#pragma once
#include "params.hpp"
#include "config.hpp"
 
#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <string>
#include <functional>
 
class DoubleQuantumDot {
public:
    /**
     * @brief Konstruktor
     * @param p parametry kropki kwantowej
     */
    explicit DoubleQuantumDot(const DotParams& p = DotParams{});
 
    /**
     * @brief Potencjał uwięzienia Vw(x) w [a.u.]
     * @param xi punkt x [a.u.]
     * @return wartość potencjału Vw(xi) [a.u.]
     * @attention ściany (|x|>a) obsługiwane przez warunki brzegowe (duże 100 na diag.)
     */
    double Vw(double xi) const;
 
    /**
     * @brief Buduje hamiltonian H bez zależności od czasu (do problemu własnego)
     * @param F_au amplituda pola elektrycznego [a.u.]
     * @return macierz hamiltonianu H(F) w bazie dyskretnej (n x n)
     */
    Eigen::MatrixXd buildH(double F_au = 0.0) const;
 
    /**
     * @brief Rozwiązuje problem własny 
     * @param F_au amplituda pola elektrycznego [a.u.]
     * @param eigvecs macierz, której kolumny to kolejne funkcje falowe (posortowanymi rosnąco wg energii)
     * @return wektor energii własnych (posortowanych rosnąco)
     */
    Eigen::VectorXd solveEigen(double F_au, Eigen::MatrixXd& eigvecs) const;
 
    /**
     * @brief Wyświetla informacje o podwójnej kropce kwantowej
     */
    inline void info() const {
        *this << std::cout; // odwrocona konwencja xD
        stream_config(std::cout);
    }
 
    friend std::ostream& operator<<(const DoubleQuantumDot& dqd, std::ostream& os) {
        os << std::fixed << dqd.get_name() <<'\n'
           << std::setprecision(0) << "  a: " << dqd.get_a_nm() << " nm, \n"
           << std::setprecision(0) << "  d1: " << dqd.get_d1_nm() << " nm, \n"
           << std::setprecision(0) << "  d2: " << dqd.get_d2_nm() << " nm, \n"
           << std::setprecision(2) << "  V1: " << dqd.get_V1_eV() << " eV, \n"
           << std::setprecision(2) << "  V2: " << dqd.get_V2_eV() << " eV, \n"
           << std::setprecision(3) << "  m_effect: " << dqd.get_m_effect() << " m_e, \n"
           << std::setprecision(0) << "  n: " << dqd.get_n() << ", \n"
           << std::setprecision(0) << "  N: " << dqd.get_N() << std::endl; 
        return os;
    }
 
    /// GETTERY 
 
    inline std::string get_name() const { return this->name_; }
    inline double get_alpha() const { return this->alpha_; }
    inline double get_dx() const { return this->dx_; }
    inline const Eigen::VectorXd& x() const { return x_; } /// zwraca wektor węzłów x [a.u.]
    
    /// FROM DOT PARAMS
    inline const DotParams& get_params() const { return this->params_; }/// Zwraca strukturę z parametrami kropki kwantowej 
    inline double get_a_nm() const { return params_.a_nm; }             /// promień kropki kwantowej w nm
    inline double get_d1_nm() const { return params_.d1_nm; }           /// zewnętrzna krawędź bariery w nm
    inline double get_d2_nm() const { return params_.d2_nm; }           /// wewnętrzna krawędź bariery w nm
    inline double get_V1_eV() const { return params_.V1_eV; }           /// głębokość bariery zewnętrznej w eV
    inline double get_V2_eV() const { return params_.V2_eV; }           /// głębokość bariery wewnętrznej w eV
    inline double get_m_effect() const { return params_.m_effect; }     /// efektywna masa w jednostkach masy elektronu
    inline int get_n() const { return params_.n; }                      /// n = liczba węzłów wewnętrznych (bez brzegów)
    inline int get_N() const { return params_.N; }                      /// N = n + 1 (z brzegami)
   

private:
    const bool fix_wavefunction_sign = true;
    static inline const std::string name_ = "Double Quantum Dot GaAs";
    DotParams params_;
    double alpha_;      /// 1/(2*m*dx^2)
    double a_au_;
    double d1_au_, d2_au_;
    double V1_au_, V2_au_;
    double dx_;
    Eigen::VectorXd x_;
 
    /**
     * @brief Inicjalizuje parametry
     */
    void initialize_parameters();
 
    /**
     * @brief Inicjalizuje siatkę x_
     */
    void initialize_grid();
 
    /**
     * @brief Normalizacja kolumn macierzy eigvecs metodą trapezów
     * @param vecs macierz do normalizacji
     * @attention Normalizacja jest konieczna, ponieważ Eigen::SelfAdjointEigenSolver nie normalizuje wektorów własnych
     */
    void norm(Eigen::MatrixXd& vecs) const;
 
};