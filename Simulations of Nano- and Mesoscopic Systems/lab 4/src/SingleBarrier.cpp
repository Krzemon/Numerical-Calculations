#include "SingleBarrier.hpp"
#include "config.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <omp.h>

SingleBarrier::SingleBarrier(double d_barrier_nm, double U_barrier_eV,
                             double m_GaAs, double m_AlGaAs, int N_slices)
    : d_barrier_nm_(d_barrier_nm)
    , U_barrier_eV_(U_barrier_eV)
    , m_GaAs_(m_GaAs)
    , m_AlGaAs_(m_AlGaAs)
    , N_slices_(N_slices)
{}

std::pair<double, double> SingleBarrier::TR_analytic(double E_eV) const
{
    const double E  = eV_to_au(E_eV);
    const double U  = eV_to_au(U_barrier_eV_);
    const double a  = nm_to_au(d_barrier_nm_);
    const double m  = m_GaAs_;  // stała masa

    double T = 0.0;

    if (std::abs(E) < tolerance || E < 0.0) {
        return {0.0, 1.0};
    }

    if (E < U) {
        // tryb tunelowy
        const double kappa = std::sqrt(2.0 * m * (U - E)) / h_bar;
        const double sh    = std::sinh(kappa * a);
        const double denom = 1.0 + (U * U * sh * sh) / (4.0 * E * (U - E));
        T = 1.0 / denom;
    } else if (std::abs(E - U) < tolerance) {
        // E ≈ U - granica
        const double ka  = std::sqrt(2.0 * m * E) / h_bar * a;
        T = 1.0 / (1.0 + m * U * U * a * a / (2.0 * h_bar * h_bar));
        (void)ka;
    } else {
        // E > U – tryb propagacyjny
        const double k  = std::sqrt(2.0 * m * (E - U)) / h_bar;
        const double sn = std::sin(k * a);
        const double denom = 1.0 + (U * U * sn * sn) / (4.0 * E * (E - U));
        T = 1.0 / denom;
    }

    return {T, 1.0 - T};
}

std::pair<double, double> SingleBarrier::TR_numeric(double E_eV, bool use_var_mass) const
{
    const double E_au = eV_to_au(E_eV);
    const double U_au = eV_to_au(U_barrier_eV_);
    const double L_contact = nm_to_au(GeometryParams::L_contact_nm);
    const double a_au = nm_to_au(d_barrier_nm_);
    const double L_au = 2.0 * L_contact + a_au;
    const int N = N_slices_;
    const double dx = L_au / static_cast<double>(N - 1);

    const double m_bar = use_var_mass ? m_AlGaAs_ : m_GaAs_;
    const cdouble k_contact = wavevec(E_au, 0.0, m_GaAs_);

    Mat2 M = mat_identity();

    for (int n = N - 2; n >= 0; --n) {
        const double x = n * dx;
        double U_struct = 0.0, m_eff = m_GaAs_;
        
        if (x < L_contact) {
            U_struct = 0.0;
            m_eff = m_GaAs_;
        } else if (x < L_contact + a_au) {
            U_struct = U_au;
            m_eff = m_bar;
        } else {
            U_struct = 0.0;
            m_eff = m_GaAs_;
        }
        
        const double x1 = (n + 1) * dx;
        double U_struct1 = 0.0, m_eff1 = m_GaAs_;
        
        if (x1 < L_contact) {
            U_struct1 = 0.0;
            m_eff1 = m_GaAs_;
        } else if (x1 < L_contact + a_au) {
            U_struct1 = U_au;
            m_eff1 = m_bar;
        } else {
            U_struct1 = 0.0;
            m_eff1 = m_GaAs_;
        }
        
        cdouble kn = wavevec(E_au, U_struct, m_eff);
        cdouble kn1 = wavevec(E_au, U_struct1, m_eff1);
        M = mat_mul(monodromy(kn, kn1, m_eff, m_eff1, x), M);
    }

    return TR_from_product(M, k_contact, k_contact, m_GaAs_, m_GaAs_);
}

std::pair<double, double> SingleBarrier::TR_const_mass(double E_eV) const
{
    return TR_numeric(E_eV, false);
}

std::pair<double, double> SingleBarrier::TR_variable_mass(double E_eV) const
{
    return TR_numeric(E_eV, true);
}

void SingleBarrier::scan_and_save(double E_min_eV, double E_max_eV, double dE_eV,
                                  const std::string& filename) const
{
    const int N = static_cast<int>((E_max_eV - E_min_eV) / dE_eV) + 1;
    std::vector<double> Es(N);
    for (int i = 0; i < N; ++i) Es[i] = E_min_eV + i * dE_eV;

    std::vector<double> Tc(N), Rc(N), Ta(N), Ra(N), Tv(N), Rv(N);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < N; ++i) {
        auto [tc, rc] = TR_const_mass(Es[i]);   Tc[i]=tc; Rc[i]=rc;
        auto [ta, ra] = TR_analytic(Es[i]);     Ta[i]=ta; Ra[i]=ra;
        auto [tv, rv] = TR_variable_mass(Es[i]); Tv[i]=tv; Rv[i]=rv;
    }

    auto fh = prepareDataFile(filename);
    stream_config(fh.out);
    fh.out << "# E[eV]  T_num_const  R_num_const  T_analytic  R_analytic  T_num_var  R_num_var\n";
    for (int i = 0; i < N; ++i)
        fh.out << Es[i] << "  " << Tc[i] << "  " << Rc[i] << "  " 
               << Ta[i] << "  " << Ra[i] << "  "
               << Tv[i] << "  " << Rv[i] << "\n";

    std::cout << " [saved]: " << fh.path << "\n";
}

void SingleBarrier::scan_var_mass_and_save(double E_min_eV, double E_max_eV, double dE_eV,
                                           const std::string& filename) const
{
    const int N = static_cast<int>((E_max_eV - E_min_eV) / dE_eV) + 1;
    std::vector<double> Es(N);
    for (int i = 0; i < N; ++i) Es[i] = E_min_eV + i * dE_eV;

    std::vector<double> Tv(N), Rv(N);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < N; ++i) {
        auto [tv, rv] = TR_variable_mass(Es[i]); Tv[i]=tv; Rv[i]=rv;
    }

    auto fh = prepareDataFile(filename);
    stream_config(fh.out);
    fh.out << "# E[eV]  T  R\n";
    for (int i = 0; i < N; ++i)
        fh.out << Es[i] << "  " << Tv[i] << "  " << Rv[i] << "\n";

    std::cout << " [saved]: " << fh.path << "\n";
}
