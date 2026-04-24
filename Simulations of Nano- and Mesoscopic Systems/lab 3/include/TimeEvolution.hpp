#pragma once
#include "params.hpp"
#include "config.hpp"
#include "DoubleQuantumDot.hpp"
 
#include <Eigen/Dense>
#include <complex>
#include <vector>
#include <ostream>
 
 
class TimeEvolution {
public:
    //  Wyniki zebrane co save_every kroków
    struct EvolutionResult {
        std::vector<double> t;      // czas 
        std::vector<double> p0;     // |psi_0|^2
        std::vector<double> p1;     // |psi_1|^2
        std::vector<double> psum;   // p0 + p1
    };
 
    using CVec = Eigen::VectorXcd; 
    using RVec = Eigen::VectorXd;
 
    /**
     * @brief Konstruktor
     * @param dot podwójna kropka kwantowa
     * @param psi0 funkcja falowa początkowa
     * @param psi1 funkcja falowa do projekcji
     * @param ep parametry ewolucji
     */
    TimeEvolution(const DoubleQuantumDot& dot,
                  const RVec& psi0, const RVec& psi1,
                  const EvolutionParams& ep = EvolutionParams{});
 
    /**
     * @brief Uruchamia ewolucję
     */
    void run();
 
    /**
     * @brief Zapisuje wyniki do strumienia wyjściowego
     */
    void write(std::ostream& os) const;
 
    /**
     * @brief Wyświetla informacje o ewolucji czasowej
     */
    inline void info() const {
        *this << std::cout; // odwrocona konwencja xD
        stream_config(std::cout);
    }
 
    friend std::ostream& operator<<(const TimeEvolution& te, std::ostream& os) {
        os << std::fixed << te.get_name() << '\n'
        << std::setprecision(0) << "  F        = " << te.get_F_kVcm()    << " [kV/cm]\n"
        << std::setprecision(0) << "  omega    = " << te.get_omega_meV() << " [meV]\n"
        << std::setprecision(0) << "  dt       = " << te.get_dt()        << " [a.u.]\n"
        << std::setprecision(0) << "  steps    = " << te.get_steps()     << "\n"
        << std::setprecision(0) << "  points   = " << te.get_result().t.size()  << "\n";
        return os;
    }

    /// GETTERY 

    inline std::string get_name() const { return this->name_; }                 /// Zwraca nazwę ewolucji czasowej
    inline const DoubleQuantumDot& get_dot() const { return this->dot_; }       /// Zwraca referencję do obiektu DoubleQuantumDot używanego w ewolucji
    inline RVec get_psi0() const { return this->psi0_; }                        /// Zwraca funkcję falową początkową 
    inline RVec get_psi1() const { return this->psi1_; }                        /// Zwraca funkcję falową do projekcji
    inline double get_F_au() const { return this->F_au_; }                      /// Zwraca amplitudę pola w jednostkach a.u.
    inline double get_omega_au() const { return this->omega_au_; }              /// Zwraca częstość w jednostkach a.u. (hbar*omega)
    inline const EvolutionResult& get_result() const { return this->result_; }  /// Zwraca wyniki ewolucji

    /// FROM EVOLUTION PARAMS
    inline const EvolutionParams& get_params() const { return params_; }
    inline double get_F_kVcm() const { return params_.F_kVcm; }
    inline double get_omega_meV() const { return params_.omega_meV; }
    inline double get_dt() const { return params_.dt; }
    inline long get_steps() const { return params_.steps; }
    inline int get_save_every() const { return params_.save_every; }
    inline int get_cn_iter() const { return params_.cn_iter; }

private:
    static inline const std::string name_ = "Quantum System Time Evolution";
    const DoubleQuantumDot& dot_;
    RVec psi0_, psi1_;
    EvolutionParams params_;
    double F_au_;
    double omega_au_;
    EvolutionResult result_;
 
    /**
     * @brief Inicjalizuje parametry
     */
    void initialize_parameters();

    /**
     * @brief Oblicza H|psi⟩ w czasie t
     * @param psi funkcja falowa
     * @param out wynik H|psi⟩
     * @param t czas
     * @attention out musi być zainicjalizowany rozmiarem N
    */
    void applyH(const CVec& psi, CVec& out, double t) const;
 
    /** 
     * @brief Krok Cranka-Nicolsona (iteracyjny)
     * @param psi_m psi^m – funkcja falowa z aktualnego kroku
     * @param t_m czas aktualnego kroku
     * @return psi^{m+1} – funkcja falowa z następnego kroku
     */
    CVec stepCN(const CVec& psi_m, double t_m) const;
 
    /**
     * @brief Krok Askara-Cakmaka (jawny)
     * @param psi_prev psi^{m-1} – funkcja falowa z poprzedniego kroku
     * @param psi_curr psi^m – funkcja falowa z aktualnego kroku
     * @param t_m czas aktualnego kroku
     * @return psi^{m+1} – funkcja falowa z następnego kroku
     */
    CVec stepAC(const CVec& psi_prev, const CVec& psi_curr, double t_m) const;
 
    /**
     * @brief Oblicza projekcję
     * @param psi funkcja falowa
     * @param ref funkcja falowa referencyjna
     * @return projekcja |<ref|psi>|^2
     */
    double proj2(const CVec& psi, const RVec& ref) const;
 
};