#include "QPCSolver.hpp"
#include "config.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <numeric>
#include <omp.h>
#include <Eigen/Eigen>

QPCSolver::QPCSolver()
    : m_eff_(QPCParams::m_eff)
    , eps_(QPCParams::eps)
    , d_nm_(QPCParams::d_nm)
    , W_nm_(QPCParams::W_nm)
    , L_nm_(QPCParams::L_nm)
    , l_frac_(QPCParams::l_frac)
    , r_frac_(QPCParams::r_frac)
    , tu_frac_(QPCParams::tu_frac)
    , bu_frac_(QPCParams::bu_frac)
    , tl_frac_(QPCParams::tl_frac)
    , bl_frac_(QPCParams::bl_frac)
    , N_subbands_(QPCParams::N_subbands)
    , Nx_(QPCParams::Nx)
    , Ny_(QPCParams::Ny)
{
    dx_ = nm_to_au(L_nm_) / static_cast<double>(Nx_ - 1);
    dy_ = nm_to_au(W_nm_) / static_cast<double>(Ny_ + 1);
    alpha_ = 0.5 / (m_eff_ * dy_ * dy_);
}

double QPCSolver::gate_potential_func(double u, double v) const
{
    const double d_au = nm_to_au(d_nm_);
    return std::atan2(u * v, d_au * std::sqrt(d_au * d_au + u * u + v * v));
}

double QPCSolver::gate_potential(double x, double y) const
{
    const double W_au = nm_to_au(W_nm_);
    const double L_au = nm_to_au(L_nm_);
    
    const double l = l_frac_ * L_au;
    const double r = r_frac_ * L_au;
    const double tu = tu_frac_ * W_au;
    const double bu = bu_frac_ * W_au;
    const double tl = tl_frac_ * W_au;
    const double bl = bl_frac_ * W_au;
    
    const double from_upper = gate_potential_func(x - l, y - bu) 
                            + gate_potential_func(x - l, tu - y)
                            + gate_potential_func(r - x, y - bu)
                            + gate_potential_func(r - x, tu - y);
    
    const double from_lower = gate_potential_func(x - l, y - bl)
                            + gate_potential_func(x - l, tl - y)
                            + gate_potential_func(r - x, y - bl)
                            + gate_potential_func(r - x, tl - y);
    
    return from_upper + from_lower;
}

std::vector<std::vector<double>> QPCSolver::compute_subbands(double V_g) const
{
    std::vector<std::vector<double>> pot_E(N_subbands_, std::vector<double>(Nx_));
    
    Eigen::VectorXd diag(Ny_);
    Eigen::VectorXd subdiag(Ny_ - 1);
    subdiag.fill(-alpha_);
    
    #pragma omp parallel for schedule(static)
    for (int nx = 0; nx < Nx_; ++nx) {
        const double x = nx * dx_;
        
        Eigen::VectorXd diag_local(Ny_);
        for (int ny = 0; ny < Ny_; ++ny) {
            const double y = ny * dy_;
            diag_local(ny) = 2.0 * alpha_ + V_g / (2.0 * M_PI * eps_) * gate_potential(x, y);
        }
        
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver;
        solver.computeFromTridiagonal(diag_local, subdiag);
        
        Eigen::VectorXd eigenvalues = solver.eigenvalues().head(N_subbands_);
        for (int n = 0; n < N_subbands_; ++n) {
            pot_E[n][nx] = eigenvalues(n);
        }
    }
    
    return pot_E;
}

void QPCSolver::save_En(const std::string& filename) const
{
    const double V_g = eV_to_au(QPCParams::V_g_eV);
    
    auto pot_E = compute_subbands(V_g);
    
    auto fh = prepareDataFile(filename);
    stream_config(fh.out);
    fh.out << "# x[nm]";
    for (int n = 0; n < N_subbands_; ++n)
        fh.out << "  E_" << n << "[meV]";
    fh.out << "\n";
    
    for (int nx = 0; nx < Nx_; ++nx) {
        fh.out << au_to_nm(nx * dx_);
        for (int n = 0; n < N_subbands_; ++n)
            fh.out << "  " << au_to_meV(pot_E[n][nx]);
        fh.out << "\n";
    }
    
    std::cout << " [saved]: " << fh.path << "\n";
}

void QPCSolver::compute_GE_and_save(double E_min_eV, double E_max_eV, double dE_eV,
                                    const std::string& filename) const
{
    const double V_g = eV_to_au(QPCParams::V_g_eV);
    
    auto pot_E = compute_subbands(V_g);
    
    const double dE_au = eV_to_au(dE_eV);
    const double E_start_au = pot_E[0][0] + dE_au;
    const double E_max_au = eV_to_au(E_max_eV);
    
    std::vector<double> Es_au;
    std::vector<double> Gs;
    
    for (double E_au = E_start_au; E_au < E_max_au; E_au += dE_au) {
        Es_au.push_back(E_au);
    }
    
    const int N_E = Es_au.size();
    Gs.resize(N_E);
    
    #pragma omp parallel for schedule(dynamic)
    for (int ie = 0; ie < N_E; ++ie) {
        const double E_au = Es_au[ie];
        
        std::vector<double> T_vec(N_subbands_);
        
        for (int n = 0; n < N_subbands_; ++n) {
            const auto& pot_En = pot_E[n];
            
            const double L_au = nm_to_au(L_nm_);
            const double dx = L_au / static_cast<double>(Nx_ - 1);
            
            const cdouble k_left = wavevec(E_au, pot_En[0], m_eff_);
            const cdouble k_right = wavevec(E_au, pot_En[Nx_ - 1], m_eff_);
            
            Mat2 M = mat_identity();
            
            for (int nx = Nx_ - 2; nx >= 0; --nx) {
                const double x = nx * dx;
                cdouble kn = wavevec(E_au, pot_En[nx], m_eff_);
                cdouble kn1 = wavevec(E_au, pot_En[nx + 1], m_eff_);
                M = mat_mul(monodromy(kn, kn1, m_eff_, m_eff_, x), M);
            }
            
            auto [T, R] = TR_from_product(M, k_left, k_right, m_eff_, m_eff_);
            T_vec[n] = T;
        }
        
        Gs[ie] = std::accumulate(T_vec.begin(), T_vec.end(), 0.0);
    }
    
    auto fh = prepareDataFile(filename);
    stream_config(fh.out);
    fh.out << "# E[meV]  G[2e^2/h]\n";
    for (int i = 0; i < N_E; ++i)
        fh.out << au_to_meV(Es_au[i]) << "  " << Gs[i] << "\n";
    
    std::cout << " [saved]: " << fh.path << "\n";
}

void QPCSolver::compute_GVg_and_save(double Vg_min_eV, double Vg_max_eV, double dVg_eV,
                                      const std::string& filename) const
{
    const double E1_au = eV_to_au(0.05);
    const double E2_au = eV_to_au(0.10);
    
    const double dVg_au = eV_to_au(dVg_eV);
    const double Vg_max_au = eV_to_au(Vg_max_eV);
    
    std::vector<double> Vgs_au;
    for (double Vg_au = 0.0; Vg_au <= Vg_max_au; Vg_au += dVg_au) {
        Vgs_au.push_back(Vg_au);
    }
    
    const int N_Vg = Vgs_au.size();
    std::vector<double> Gs1(N_Vg);
    std::vector<double> Gs2(N_Vg);
    
    for (int ivg = 0; ivg < N_Vg; ++ivg) {
        const double V_g = Vgs_au[ivg];
        
        auto pot_E = compute_subbands(V_g);
        
        std::vector<double> T_vec1(N_subbands_);
        #pragma omp parallel for
        for (int n = 0; n < N_subbands_; ++n) {
            const auto& pot_En = pot_E[n];
            
            const double L_au = nm_to_au(L_nm_);
            const double dx = L_au / static_cast<double>(Nx_ - 1);
            const cdouble k_left = wavevec(E1_au, pot_En[0], m_eff_);
            const cdouble k_right = wavevec(E1_au, pot_En[Nx_ - 1], m_eff_);
            
            Mat2 M = mat_identity();
            for (int nx = Nx_ - 2; nx >= 0; --nx) {
                const double x = nx * dx;
                cdouble kn = wavevec(E1_au, pot_En[nx], m_eff_);
                cdouble kn1 = wavevec(E1_au, pot_En[nx + 1], m_eff_);
                M = mat_mul(monodromy(kn, kn1, m_eff_, m_eff_, x), M);
            }
            
            auto [T, R] = TR_from_product(M, k_left, k_right, m_eff_, m_eff_);
            T_vec1[n] = T;
        }
        Gs1[ivg] = std::accumulate(T_vec1.begin(), T_vec1.end(), 0.0);
        
        std::vector<double> T_vec2(N_subbands_);
        #pragma omp parallel for
        for (int n = 0; n < N_subbands_; ++n) {
            const auto& pot_En = pot_E[n];
            
            const double L_au = nm_to_au(L_nm_);
            const double dx = L_au / static_cast<double>(Nx_ - 1);
            const cdouble k_left = wavevec(E2_au, pot_En[0], m_eff_);
            const cdouble k_right = wavevec(E2_au, pot_En[Nx_ - 1], m_eff_);
            
            Mat2 M = mat_identity();
            for (int nx = Nx_ - 2; nx >= 0; --nx) {
                const double x = nx * dx;
                cdouble kn = wavevec(E2_au, pot_En[nx], m_eff_);
                cdouble kn1 = wavevec(E2_au, pot_En[nx + 1], m_eff_);
                M = mat_mul(monodromy(kn, kn1, m_eff_, m_eff_, x), M);
            }
            
            auto [T, R] = TR_from_product(M, k_left, k_right, m_eff_, m_eff_);
            T_vec2[n] = T;
        }
        Gs2[ivg] = std::accumulate(T_vec2.begin(), T_vec2.end(), 0.0);
    }
    
    auto fh = prepareDataFile(filename);
    stream_config(fh.out);
    fh.out << "# V_g[meV]  G_50meV[2e^2/h]  G_100meV[2e^2/h]\n";
    for (int i = 0; i < N_Vg; ++i)
        fh.out << au_to_meV(Vgs_au[i]) << "  " << Gs1[i] << "  " << Gs2[i] << "\n";
    
    std::cout << " [saved]: " << fh.path << "\n";
}
