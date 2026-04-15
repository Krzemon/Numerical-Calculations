#include "ImaginaryTime.hpp"

#include <iomanip>
#include <cmath>
#include <random>
#include <iostream>


ImaginaryTime::ImaginaryTime(double a_nm_in, int N_in, double tol_in): 
    N_(N_in), tol_(tol_in)
{
    initialize_parameters(a_nm_in);
    initialize_grid();
    coulomb_interaction();
}


std::ostream& operator<<(std::ostream& os, const ImaginaryTime& it)
{
    os << "  Zbieżność osiągnięta po "
       << it.iter_ << " iteracjach\n";

    os << "  Energia stanu podstawowego: "
       << it.E_final_ << " [eV]\n";

    os << "  Liczba węzłów: "
       << it.N_ << "\n";

    os << "  a: "
       << std::fixed << std::setprecision(1)
       << it.a_au_ * bohr_ratio_nm
       << " [nm]\n\n";

    stream_config();

    return os;
}


void ImaginaryTime::run()
{
    compute();
}


void ImaginaryTime::write_energies(std::ostream& os) const
{
    for (size_t k = 0; k < E_history_.size(); ++k)
        os << k*100 << '\t' << E_history_[k] << '\n';
}


void ImaginaryTime::write_density(std::ostream& os, bool is_nanometers) const
{
    os << "# units=" << (is_nanometers ? "nm" : "au") << "\n";

    for (int i = 0; i < N_; ++i)
    {
        double x1 = -a_au_ + i * dx_;
        if (is_nanometers) x1 *= bohr_ratio_nm;

        for (int j = 0; j < N_; ++j)
        {
            double x2   = -a_au_ + j * dx_;
            double psi2 = Psi_[index(i,j)] * Psi_[index(i,j)];

            if (is_nanometers)
            {
                x2   *= bohr_ratio_nm;
                psi2 /= (bohr_ratio_nm * bohr_ratio_nm);
            }

            os << x1 << '\t' << x2 << '\t' << psi2 << '\n';
        }

        os << '\n';
    }
}

/// PRIVATE

void ImaginaryTime::initialize_parameters(double a_nm_in)
{
    a_au_ = a_nm_in / bohr_ratio_nm;
    l_au_ = l_nm / bohr_ratio_nm;
    dx_ = 2.0 * a_au_ / (N_ - 1);
    dt_ = 0.4 * m_effect * dx_ * dx_;
}


void ImaginaryTime::initialize_grid()
{
    Psi_.assign(N_*N_, 0.0);
    V_int_.resize(N_);
}


void ImaginaryTime::coulomb_interaction()
{
    for (int k = 0; k < N_; ++k)
    {
        double r = k * dx_;
        V_int_[k] = 1.0 / (eps_r * std::sqrt(r*r + l_au_*l_au_));
    }
}


void ImaginaryTime::initialize_random()
{
    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<double> dist(-1.0, 1.0);

    for (int i = 0; i < N_; ++i)
    {
        for (int j = 0; j < N_; ++j)
        {
            if (is_boundary(i,j))
                Psi_[index(i,j)] = 0.0;
            else
                Psi_[index(i,j)] = dist(rng);
        }
    }
    normalize(Psi_);
}


void ImaginaryTime::apply_hamiltonian(const Grid& psi, Grid& Hpsi) const
{
    const double alpha = 1.0 / (2.0 * m_effect * dx_ * dx_);
    std::fill(Hpsi.begin(), Hpsi.end(), 0.0);

    for (int i = 1; i < N_-1; ++i)
    {
        for (int j = 1; j < N_-1; ++j)
        {
            double laplace =
                    psi[index(i+1,j)]
                + psi[index(i-1,j)]
                + psi[index(i,j+1)]
                + psi[index(i,j-1)]
                - 4.0 * psi[index(i,j)];

            double kinetic = -alpha * laplace;
            double potential = V_int_[std::abs(i-j)] * psi[index(i,j)];
            Hpsi[index(i,j)] = kinetic + potential;
        }
    }
}


void ImaginaryTime::normalize(Grid& psi) const
{
    double norm2 = 0.0;
    for (int i = 1; i < N_-1; ++i)
    {
        for (int j = 1; j < N_-1; ++j)
        {
            double v = psi[index(i,j)];
            norm2 += v*v;
        }
    }

    norm2 *= dx_*dx_;
    double inv = 1.0 / std::sqrt(norm2);
    for (double& v : psi)
        v *= inv;
}


double ImaginaryTime::compute_energy(const Grid& psi, const Grid& Hpsi) const
{
    double E = 0.0;
    for (int i = 1; i < N_-1; ++i)
        for (int j = 1; j < N_-1; ++j)
            E += psi[index(i,j)] * Hpsi[index(i,j)];

    return E * dx_ * dx_;
}


void ImaginaryTime::compute()
{
    initialize_random();
    Grid Hpsi(N_*N_, 0.0);

    double E = 0.0;
    double E_prev = 1.0;
    long long iter = 0;

    while (std::abs((E - E_prev)/E_prev) > tol_)
    {
        ++iter;
        E_prev = E;

        apply_hamiltonian(Psi_, Hpsi); // |H|psi>

        // Imaginary-time step
        for (int i = 1; i < N_-1; ++i)
            for (int j = 1; j < N_-1; ++j)
                Psi_[index(i,j)] -= dt_ * Hpsi[index(i,j)];

        normalize(Psi_);
        apply_hamiltonian(Psi_, Hpsi);
        E = compute_energy(Psi_, Hpsi);

        if (iter % 100 == 0)
        {
            double E_eV = E * E_h;
            E_history_.push_back(E_eV);
        }
    }
    E_final_ = E*E_h;
    iter_ = iter;
}