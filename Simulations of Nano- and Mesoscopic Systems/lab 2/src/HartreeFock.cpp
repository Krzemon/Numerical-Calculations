#include "HartreeFock.hpp"

#include <iomanip>
#include <cmath>
#include <random>
#include <iostream>
 

HartreeFock::HartreeFock(double a_nm_in, int N_in,double tol_IT_in, double tol_HF_in): 
N_(N_in), tol_IT_(tol_IT_in), tol_HF_(tol_HF_in)
{
    initialize_parameters(a_nm_in);
    initialize_grid();
    coulomb_interaction();
}


std::ostream& operator<<(std::ostream& os, const HartreeFock& hf)
{
    os << "  Zbieżność osiągnięta po "
       << hf.iter_ << " iteracjach\n";

    os << "  Energia stanu podstawowego: "
       << hf.E_final_ << " [eV]\n";

    os << "  Liczba węzłów: "
       << hf.N_ << "\n";

    os << "  a: "
       << std::fixed << std::setprecision(1)
       << hf.a_au_ * bohr_ratio_nm
       << " [nm]\n\n";

    stream_config();

    return os;
}


void HartreeFock::run()
{
    compute();
}


void HartreeFock::write_energies(std::ostream& os) const
{
    for (size_t k = 0; k < E_history_.size(); ++k)
        os << k << '\t' << E_history_[k] << '\n';
}


void HartreeFock::write_density(std::ostream& os, bool is_nanometers) const
{
    os << "# units=" << (is_nanometers ? "nm" : "au") << "\n";

    for (int i = 0; i < N_; ++i) {
        double x    = -a_au_ + i * dx_;
        double psi2 = Psi_[i] * Psi_[i];
        if (is_nanometers) {
            x *= bohr_ratio_nm;
            psi2 /= bohr_ratio_nm;
        }
        os << x << '\t' << psi2 << '\n';
    }
}

/// PRIVATE

void HartreeFock::initialize_parameters(double a_nm_in)
{
    a_au_ = a_nm_in / bohr_ratio_nm;
    l_au_ = l_nm / bohr_ratio_nm;
    dx_ = 2.0 * a_au_ / (N_ - 1);
    dt_ = 0.4 * m_effect * dx_ * dx_;
}


void HartreeFock::initialize_grid()
{
    Psi_.assign(N_, 0.0);
    J1_.assign(N_, 0.0);
    V_int_.resize(N_);
}


void HartreeFock::coulomb_interaction()
{
    for (int k = 0; k < N_; ++k) {
        double r = k * dx_;
        V_int_[k] = 1.0 / (eps_r * std::sqrt(r*r + l_au_*l_au_));
    }
}


void HartreeFock::initialize_random()
{ 
    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    for (int i = 0; i < N_; ++i)
        Psi_[i] = is_boundary(i) ? 0.0 : dist(rng);
    normalize(Psi_);
}


void HartreeFock::normalize(Grid& psi) const
{
    double norm2 = 0.0;
    for (int i = 1; i < N_-1; ++i)
        norm2 += psi[i] * psi[i];
    norm2 *= dx_;
    double inv = 1.0 / std::sqrt(norm2);
    for (double& v : psi) 
        v *= inv;
}


void HartreeFock::compute_J1()
{
    std::fill(J1_.begin(), J1_.end(), 0.0);
    for (int i = 1; i < N_-1; ++i) {
        double sum = 0.0;
        for (int j = 1; j < N_-1; ++j)
            sum += Psi_[j] * Psi_[j] * V_int_[std::abs(i-j)];
        J1_[i] = sum * dx_;
    }
}


void HartreeFock::apply_fock(const Grid& psi, Grid& Fpsi) const
{
    const double alpha = 1.0 / (2.0 * m_effect * dx_ * dx_);
    std::fill(Fpsi.begin(), Fpsi.end(), 0.0);
    for (int i = 1; i < N_-1; ++i) {
        double kinetic = -alpha * (psi[i+1] + psi[i-1] - 2.0*psi[i]);
        double hartree = J1_[i] * psi[i];
        Fpsi[i] = kinetic + hartree;
    }
}


double HartreeFock::compute_epsilon(const Grid& psi, const Grid& Fpsi) const
{
    double eps = 0.0;
    for (int i = 1; i < N_-1; ++i)
        eps += psi[i] * Fpsi[i];
    return eps * dx_;
}


double HartreeFock::compute_E_int() const
{
    double E = 0.0;
    for (int i = 1; i < N_-1; ++i)
        E += Psi_[i] * Psi_[i] * J1_[i];
    return E * dx_;
}


double HartreeFock::solve_fock_imaginary_time()
{
    Grid Fpsi(N_, 0.0);

    double eps = 0.0;
    double eps_old = 1.0;

    while (std::abs((eps - eps_old) / (std::abs(eps_old) + tol_IT_)) > tol_IT_)
    {
        eps_old = eps;

        apply_fock(Psi_, Fpsi);

        for (int i = 1; i < N_-1; ++i)
            Psi_[i] -= dt_ * Fpsi[i];

        normalize(Psi_);

        apply_fock(Psi_, Fpsi);
        eps = compute_epsilon(Psi_, Fpsi);
    }

    return eps;
}


void HartreeFock::compute()
{
    initialize_random();

    double E = 0.0;
    double E_prev = 1.0;
    long long iter = 0;

    while (std::abs((E - E_prev) / (std::abs(E_prev) + tol_HF_)) > tol_HF_)
    {
        ++iter;
        E_prev = E;

        compute_J1();
        double eps1 = solve_fock_imaginary_time();
        double E_int = compute_E_int();
        E = 2.0 * eps1 - E_int;

        E_history_.push_back(E * E_h);
    }
    E_final_ = E * E_h;
    iter_ = iter;
}