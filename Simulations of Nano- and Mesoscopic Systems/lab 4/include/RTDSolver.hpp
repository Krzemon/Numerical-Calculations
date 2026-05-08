#pragma once

#include "TransferMatrix.hpp"
#include "params.hpp"

#include <vector>
#include <string>

/**
 * @brief Rozwiązuje transport przez diodę rezonansowo-tunelowaną (RTD).
 */
class RTDSolver : public TransferMatrix {
public:
    explicit RTDSolver(double d_barrier_nm = RTDParams::d_barrier_nm,
                       double d_well_nm = RTDParams::d_well_nm,
                       double U_barrier_eV = RTDParams::U_barrier_eV,
                       double m_GaAs = RTDParams::m_GaAs,
                       double m_AlGaAs = calc_m_AlGaAs(RTDParams::m_AlGaAs_x),
                       int N_slices = RTDParams::N_slices);

    std::pair<double, double> TR(double E_eV, double Vbias = 0.0) const;
    
    /**
     * @brief Skanuje współczynnik transmisji i odbicia dla zakresu energii, a następnie zapisuje wyniki do pliku.
     */
    void scan_TE_and_save(double E_min_eV, double E_max_eV, double dE_eV,
                          const std::string& filename = "zad3_rtd_TE.dat") const;
    /**
     * @brief Oblicza charakterystykę prądowo-napięciową, a następnie zapisuje wyniki do pliku.
     */
    void compute_IV_and_save(double Vbias_min = 0.0, double Vbias_max = 0.5, double dV = 0.001,
                             double mu_eV = RTDParams::mu_eV, const std::string& filename = "zad3_rtd_IV.dat") const;

    /**
     * @brief Zwraca całkowitą długość struktury RTD w nanometrach.
     */
    double total_length_nm() const;

private:
    double d_barrier_nm_;
    double d_well_nm_;
    double U_barrier_eV_;
    double m_GaAs_;
    double m_AlGaAs_;
    int N_slices_;

    std::pair<double, double> TR_numeric(double E_eV, double Vbias) const;
};
