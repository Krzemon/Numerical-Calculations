#pragma once

#include "TransferMatrix.hpp"
#include "params.hpp"

#include <vector>
#include <string>
#include <functional>

/**
 * @brief Rozwiązuje transport przez kwantowy kontakt punktowy (QPC).
 */
class QPCSolver : public TransferMatrix {
public:
    explicit QPCSolver();

    /**
     * @brief Oblicza i zapisuje profile energii do pliku.
     */
    void save_En(const std::string& filename = "zad4_En.dat") const;
    
    /**
     * @brief Oblicza przewodność GE dla zakresu energii, a następnie zapisuje wyniki do plików.
     */
    void compute_GE_and_save(double E_min_eV, double E_max_eV, double dE_eV, 
                             const std::string& filename = "zad4_GE.dat") const;
    
    /**
     * @brief Oblicza przewodność GVg dla zakresu potencjałów bramki, a następnie zapisuje wyniki do pliku.
     */
    void compute_GVg_and_save(double Vg_min_eV, double Vg_max_eV, double dVg_eV, 
                              const std::string& filename = "zad4_GVg.dat") const;

private:
    double m_eff_;
    double eps_;
    double d_nm_;
    double W_nm_;
    double L_nm_;
    
    double l_frac_, r_frac_;
    double tu_frac_, bu_frac_;
    double tl_frac_, bl_frac_;
    
    int N_subbands_;
    int Nx_;
    int Ny_;
    
    double dx_;
    double dy_;
    double alpha_;
    
    double gate_potential_func(double u, double v) const;
    double gate_potential(double x, double y) const;
    std::vector<std::vector<double>> compute_subbands(double V_g) const;
};
