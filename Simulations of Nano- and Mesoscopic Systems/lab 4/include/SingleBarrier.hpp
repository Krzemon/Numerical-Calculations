#pragma once

#include "TransferMatrix.hpp"
#include "params.hpp"

#include <vector>
#include <string>

/**
 * @brief Klasa rozwiązująca transport przez pojedynczą barierę 
 *        GaAs/AlGaAs metodą macierzy transferu.
 */
class SingleBarrier : public TransferMatrix {
public:
    /**
     * @brief Konstruktor – przyjmuje parametry urządzenia.
     */
    explicit SingleBarrier(double d_barrier_nm  = 5.0, double U_barrier_eV = 0.27, double m_GaAs = 0.063,
                           double m_AlGaAs = calc_m_AlGaAs(0.3), int N_slices = 500);

    /**
     * @brief Współczynnik transmisji T(E) przy stałej masie (Zadanie 1).
     */
    std::pair<double, double> TR_const_mass(double E_eV) const;

    /**
     * @brief Współczynnik transmisji T(E) przy zmiennej masie (Zadanie 2).
     */
    std::pair<double, double> TR_variable_mass(double E_eV) const;

    /**
     * @brief Analityczny współczynnik transmisji dla prostokątnej bariery.
     */
    std::pair<double, double> TR_analytic(double E_eV) const;

    /**
     * @brief Skanuje zakres energii i zapisuje wyniki do pliku (zadanie 1).
     */
    void scan_and_save(double E_min_eV, double E_max_eV, double dE_eV,
                       const std::string& filename = "zad1_single_barrier.dat") const;

    /**
     * @brief Skanuje zakres energii dla zmiennej masy (zadanie 2).
     */
    void scan_var_mass_and_save(double E_min_eV, double E_max_eV, double dE_eV,
                                const std::string& filename = "zad2_single_barrier_var_mass.dat") const;

    // gettery
    double get_d_nm() const { return d_barrier_nm_; }
    double get_U_eV() const { return U_barrier_eV_; }
    double get_m_GaAs() const { return m_GaAs_; }
    double get_m_AlGaAs() const { return m_AlGaAs_; }

private:
    double d_barrier_nm_;
    double U_barrier_eV_;
    double m_GaAs_;
    double m_AlGaAs_;
    int N_slices_;

    /**
     * @brief Wspólna logika numeryczna – buduje profil i mnoży macierze.
     */
    std::pair<double, double> TR_numeric(double E_eV, bool use_var_mass) const;
};
