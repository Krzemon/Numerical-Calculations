#include "TimeEvolution.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
 
using namespace std::complex_literals; // liczby zespolone w notacji a + bi
 

TimeEvolution::TimeEvolution(const DoubleQuantumDot& dot, 
                             const RVec& psi0, 
                             const RVec& psi1, 
                             const EvolutionParams& ep)
    : dot_(dot), psi0_(psi0), psi1_(psi1), params_(ep)
{
    initialize_parameters();
}


void TimeEvolution::initialize_parameters()
{
    this->F_au_ = F_kVcm_to_au(get_F_kVcm());
    this->omega_au_ = meV_to_au(get_omega_meV());
}


void TimeEvolution::applyH(const CVec& psi, CVec& out, double t) const
{
    const double alpha = dot_.get_alpha();
    const double Ft = get_F_au() * std::sin(get_omega_au() * t);
    const auto& x = dot_.x();
    const int n = dot_.get_params().n;
 
    for (int i = 0; i < n; ++i) {
        std::complex<double> left  = (i > 0)   ? psi(i-1) : 0.0;
        std::complex<double> right = (i < n-1) ? psi(i+1) : 0.0;
        out(i) = -alpha*(right + left - 2.0*psi(i)) + dot_.Vw(x(i))*psi(i) + Ft*x(i)*psi(i);
    }
}
 
 
TimeEvolution::CVec TimeEvolution::stepCN(const CVec& psi_m, double t_m) const
{
    const std::complex<double> coeff = params_.dt / (2.0 * 1i);
    const double t_m1 = t_m + params_.dt;
    const int n = dot_.get_params().n;
 
    CVec Hm_psim(n);
    applyH(psi_m, Hm_psim, t_m);
    CVec psi_k = psi_m;
    CVec psi_k_prev(n);
 
    do {
        psi_k_prev = psi_k;
        CVec Hm1_psik(n);
        applyH(psi_k, Hm1_psik, t_m1);
        psi_k = psi_m + coeff * (Hm_psim + Hm1_psik);
    } while ((psi_k - psi_k_prev).norm() > tolerance);
 
    return psi_k;
}
 
 
TimeEvolution::CVec TimeEvolution::stepAC(const CVec& psi_prev, const CVec& psi_curr, double t_m) const
{
    // Askar-Cakmak: psi_{m+1} = psi_{m-1} - 2i*dt * H(t_m) * psi_m, t_m = m*dt
    const std::complex<double> coeff = 2.0 * params_.dt / 1i;
    CVec Hpsi(dot_.get_params().n);
    applyH(psi_curr, Hpsi, t_m);
    return psi_prev + coeff * Hpsi;
}
 
 
double TimeEvolution::proj2(const CVec& psi, const RVec& ref) const
{
    // <ref|psi> = sum_i ref(i)*psi(i)*dx
    // |<ref|psi>|^2 = norm(overlap)^2 = norm(overlap * dx)^2 = norm(overlap)^2 * dx^2
    const double dx = dot_.get_dx();
    std::complex<double> overlap = ref.cast<std::complex<double>>().dot(psi);
    return std::norm(overlap) * dx * dx;
}
 
 
void TimeEvolution::run()
{
    result_ = {};   // czyść poprzednie wyniki
 
    const int nout = params_.steps / params_.save_every + 2;
    result_.t.reserve(nout);
    result_.p0.reserve(nout);
    result_.p1.reserve(nout);
    result_.psum.reserve(nout);
 
    // Zapis t=0
    CVec psi_curr = psi0_.cast<std::complex<double>>();
    auto save = [&](double t, const CVec& psi) {
        double p0 = proj2(psi, psi0_);
        double p1 = proj2(psi, psi1_);
        result_.t.push_back(t);
        result_.p0.push_back(p0);
        result_.p1.push_back(p1);
        result_.psum.push_back(p0 + p1);
    };
    save(0.0, psi_curr);
 
    // Crank-Nicolson (pierwszy krok)
    CVec psi_prev = psi_curr;
    psi_curr = stepCN(psi_curr, 0.0);
 
    // Askar-Cakmak
    for (int m = 1; m < params_.steps; ++m) {
        if (m % params_.save_every == 0)
            save(m * params_.dt, psi_curr);
 
        const double t_m = m * params_.dt;
        CVec psi_next = stepAC(psi_prev, psi_curr, t_m);
        psi_prev = psi_curr;
        psi_curr = psi_next;
 
        if (m % 100'000 == 0) {
            double p0 = proj2(psi_curr, psi0_);
            double p1 = proj2(psi_curr, psi1_);
            std::cerr << "  [AC] step = " << std::setw(8) << m
                      << "  t = " << std::setw(10) << au_to_ns(m * params_.dt) << " ns"
                      << "  p0 = " << std::setw(8) << std::setprecision(4) << p0
                      << "  p1 = " << std::setw(8) << p1
                      << "  sum = " << (p0 + p1) << "\n";
        }
    }
}
 
 
void TimeEvolution::write(std::ostream& os) const
{
    stream_config(os);
    os << "# t[ns] |<Psi|0>|^2  |<Psi|1>|^2  sum\n";
    for (std::size_t i = 0; i < result_.t.size(); ++i)
        os << au_to_ns(result_.t[i]) << "  " << result_.p0[i] << "  "
           << result_.p1[i]<< "  " << result_.psum[i] << "\n";
}
 
 
std::ostream& operator<<(std::ostream& os, const TimeEvolution& te)
{
    const auto& ep = te.get_params();
    const auto& r  = te.get_result();
    os << "TimeEvolution statystyki:\n"
       << "  F        = " << ep.F_kVcm   << " [kV/cm]\n"
       << "  omega    = " << ep.omega_meV << " [meV]\n"
       << "  dt       = " << ep.dt       << " [a.u.]\n"
       << "  steps    = " << ep.steps   << "\n"
       << "  points   = " << r.t.size()  << "\n";
    if (!r.psum.empty())
        os << "  psum(end)= " << r.psum.back()
           << "  (oczekiwane ≈ 1)\n";
    return os;
}