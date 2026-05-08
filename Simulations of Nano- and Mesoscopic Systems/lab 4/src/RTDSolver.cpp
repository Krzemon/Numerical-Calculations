#include "RTDSolver.hpp"
#include "config.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <numeric>
#include <omp.h>

RTDSolver::RTDSolver(double d_barrier_nm, double d_well_nm, double U_barrier_eV, 
                     double m_GaAs, double m_AlGaAs, int N_slices)
    : d_barrier_nm_(d_barrier_nm)
    , d_well_nm_(d_well_nm)
    , U_barrier_eV_(U_barrier_eV)
    , m_GaAs_(m_GaAs)
    , m_AlGaAs_(m_AlGaAs)
    , N_slices_(N_slices)
{}

double RTDSolver::total_length_nm() const
{
    return 2.0 * GeometryParams::L_contact_nm + 2.0 * d_barrier_nm_ + d_well_nm_;
}

std::pair<double, double> RTDSolver::TR_numeric(double E_eV, double Vbias) const
{
    const double E_au = eV_to_au(E_eV);
    const double U_au = eV_to_au(U_barrier_eV_);
    const double Vbias_au = eV_to_au(Vbias);
    
    const double L_contact = nm_to_au(GeometryParams::L_contact_nm);
    const double d_bar_au = nm_to_au(d_barrier_nm_);
    const double d_well_au = nm_to_au(d_well_nm_);
    const double L_au = 2.0 * L_contact + 2.0 * d_bar_au + d_well_au;
    
    const int N = N_slices_;
    const double dx = L_au / static_cast<double>(N - 1);

    const cdouble k_left = wavevec(E_au, 0.0, m_GaAs_);
    const cdouble k_right = wavevec(E_au, -Vbias_au, m_GaAs_);

    Mat2 M = mat_identity();

    for (int n = N - 2; n >= 0; --n) {
        const double x = n * dx;
        const double x1 = (n + 1) * dx;
        
        const double U_bias = -Vbias_au * x / L_au;
        const double U_bias1 = -Vbias_au * x1 / L_au;
        
        double U_struct = 0.0, m_eff = m_GaAs_;
        if (x < L_contact) {
            U_struct = 0.0;
            m_eff = m_GaAs_;
        } else if (x < L_contact + d_bar_au) {
            U_struct = U_au;
            m_eff = m_AlGaAs_;
        } else if (x < L_contact + d_bar_au + d_well_au) {
            U_struct = 0.0;
            m_eff = m_GaAs_;
        } else if (x < L_contact + 2.0 * d_bar_au + d_well_au) {
            U_struct = U_au;
            m_eff = m_AlGaAs_;
        } else {
            U_struct = 0.0;
            m_eff = m_GaAs_;
        }
        
        double U_struct1 = 0.0, m_eff1 = m_GaAs_;
        if (x1 < L_contact) {
            U_struct1 = 0.0;
            m_eff1 = m_GaAs_;
        } else if (x1 < L_contact + d_bar_au) {
            U_struct1 = U_au;
            m_eff1 = m_AlGaAs_;
        } else if (x1 < L_contact + d_bar_au + d_well_au) {
            U_struct1 = 0.0;
            m_eff1 = m_GaAs_;
        } else if (x1 < L_contact + 2.0 * d_bar_au + d_well_au) {
            U_struct1 = U_au;
            m_eff1 = m_AlGaAs_;
        } else {
            U_struct1 = 0.0;
            m_eff1 = m_GaAs_;
        }
        
        const double U_total = U_struct + U_bias;
        const double U_total1 = U_struct1 + U_bias1;
        
        cdouble kn = wavevec(E_au, U_total, m_eff);
        cdouble kn1 = wavevec(E_au, U_total1, m_eff1);
        M = mat_mul(monodromy(kn, kn1, m_eff, m_eff1, x), M);
    }

    return TR_from_product(M, k_left, k_right, m_GaAs_, m_GaAs_);
}

std::pair<double, double> RTDSolver::TR(double E_eV, double Vbias) const
{
    return TR_numeric(E_eV, Vbias);
}

void RTDSolver::scan_TE_and_save(double E_min_eV, double E_max_eV, double dE_eV,
                                 const std::string& filename) const
{
    const int N = static_cast<int>((E_max_eV - E_min_eV) / dE_eV) + 1;
    std::vector<double> Es(N);
    for (int i = 0; i < N; ++i) Es[i] = E_min_eV + i * dE_eV;

    std::vector<double> Ts(N), Rs(N);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < N; ++i) {
        auto [T, R] = TR(Es[i], 0.0);
        Ts[i] = T;
        Rs[i] = R;
    }

    auto fh = prepareDataFile(filename);
    stream_config(fh.out);
    fh.out << "# E[eV]  T  R\n";
    for (int i = 0; i < N; ++i)
        fh.out << Es[i] << "  " << Ts[i] << "  " << Rs[i] << "\n";

    std::cout << " [saved]: " << fh.path << "\n";
}

void RTDSolver::compute_IV_and_save(double Vbias_min, double Vbias_max, double dV, double mu_eV, 
                                    const std::string& filename) const
{
    const double kB_eV_K = 8.617333262e-5;
    const double T_He_K = RTDParams::T_He_K;
    const double T_N_K = RTDParams::T_N_K;
    
    const double kBT_He_eV = kB_eV_K * T_He_K;
    const double kBT_N_eV = kB_eV_K * T_N_K;
    
    const double kBT_He_au = eV_to_au(kBT_He_eV);
    const double kBT_N_au = eV_to_au(kBT_N_eV);
    
    const double mu_s_au = eV_to_au(mu_eV);
    const double mu_d_au = eV_to_au(mu_eV);
    
    const int N_V = static_cast<int>((Vbias_max - Vbias_min) / dV) + 1;
    std::vector<double> Vs(N_V);
    for (int i = 0; i < N_V; ++i) Vs[i] = Vbias_min + i * dV;

    std::vector<double> js_He(N_V);
    std::vector<double> js_N(N_V);

    #pragma omp parallel for schedule(dynamic)
    for (int iv = 0; iv < N_V; ++iv) {
        const double Vbias = Vs[iv];
        const double Vbias_au = eV_to_au(Vbias);
        
        std::vector<double> jE_He_vec;
        std::vector<double> jE_N_vec;
        
        const double dE_au = eV_to_au(0.0001);
        const double E_max_au = eV_to_au(0.2);
        const int N_E = static_cast<int>(E_max_au / dE_au);
        jE_He_vec.reserve(N_E + 1);
        jE_N_vec.reserve(N_E + 1);
        
        for (double E_au = dE_au; E_au < E_max_au; E_au += dE_au) {
            const double E_eV = au_to_eV(E_au);
            
            auto [T, R] = TR(E_eV, Vbias);
            
            const double log_term_He = std::log(
                (1.0 + std::exp((mu_s_au - E_au) / kBT_He_au)) / 
                (1.0 + std::exp((mu_d_au - Vbias_au - E_au) / kBT_He_au))
            );
            
            const double log_term_N = std::log(
                (1.0 + std::exp((mu_s_au - E_au) / kBT_N_au)) / 
                (1.0 + std::exp((mu_d_au - Vbias_au - E_au) / kBT_N_au))
            );
            
            jE_He_vec.push_back(1e6 * T * log_term_He);
            jE_N_vec.push_back(1e6 * T * log_term_N);
        }
        
        js_He[iv] = std::accumulate(jE_He_vec.begin(), jE_He_vec.end(), 0.0) * dE_au * kBT_He_au;
        js_N[iv] = std::accumulate(jE_N_vec.begin(), jE_N_vec.end(), 0.0) * dE_au * kBT_N_au;
    }

    auto fh = prepareDataFile(filename);
    stream_config(fh.out);
    fh.out << "# Vbias[mV]  j_He[a.u.]  j_N[a.u.]\n";
    for (int i = 0; i < N_V; ++i)
        fh.out << Vs[i] * 1000.0 << "  " << js_He[i] << "  " << js_N[i] << "\n";

    std::cout << " [saved]: " << fh.path << "\n";
}
