#include "DoubleQuantumDot.hpp"
#include <Eigen/Eigenvalues>
#include <cmath>
#include <algorithm>
#include <numeric>
 
DoubleQuantumDot::DoubleQuantumDot(const DotParams& p)
    : params_(p)
{
    initialize_parameters();
    initialize_grid();
}
 

double DoubleQuantumDot::Vw(double xi) const
{
    double ax = std::abs(xi);
    if (ax <= d2_au_)            // bariera wewnętrzna
        return V2_au_;
    if (ax <= a_au_ && ax >= d1_au_) // bariera zewnętrzna
        return V1_au_;
    return 0.0;
}


Eigen::MatrixXd DoubleQuantumDot::buildH(double F_au) const
{
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(get_n(), get_n());
 
    for (int i = 0; i < get_n(); ++i) {
        H(i, i) = 2.0 * alpha_ + Vw(x_(i)) + F_au * x_(i);
        if (i > 0)         H(i, i-1) = -alpha_;
        if (i < get_n()-1) H(i, i+1) = -alpha_;
    }
    return H;
}


void DoubleQuantumDot::initialize_parameters()
{
    this->a_au_  = nm_to_au(get_a_nm());
    this->d1_au_ = nm_to_au(get_d1_nm());
    this->d2_au_ = nm_to_au(get_d2_nm());
    this->V1_au_ = eV_to_au(get_V1_eV());
    this->V2_au_ = eV_to_au(get_V2_eV());
    this->dx_ = 2.0 * a_au_ / static_cast<double>(get_n() + 1);
    this->alpha_ = 1.0 / (2.0 * get_m_effect() * dx_*dx_);
}


void DoubleQuantumDot::initialize_grid()
{
    this->x_.resize(get_n());
    for (int i = 0; i < get_n(); ++i)
        this->x_(i) = -a_au_ + (i + 1) * dx_;
}


void DoubleQuantumDot::norm(Eigen::MatrixXd& vecs) const
{
    for (int k = 0; k < vecs.cols(); ++k) 
    {
        double norm2 = 0.0;
        for (int i = 0; i < get_n(); ++i)
            norm2 += vecs(i, k) * vecs(i, k);
        norm2 *= dx_;
        double norm = std::sqrt(norm2);
        if (norm > 1e-14)
            vecs.col(k) /= norm;
    }
}
 

Eigen::VectorXd DoubleQuantumDot::solveEigen(double F_au, Eigen::MatrixXd& eigvecs) const
{
    Eigen::MatrixXd H = buildH(F_au);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(H);
    if (solver.info() != Eigen::Success)
        throw std::runtime_error("Eigen: diagonalizacja nie powiodla sie!");
 
    eigvecs = solver.eigenvectors();  // kolumny posortowane rosnąco wg E
    norm(eigvecs);
 
    if(fix_wavefunction_sign) {
        for (int k = 0; k < eigvecs.cols(); ++k) {
            int imax = 0;
            double vmax = 0.0;
            for (int i = 0; i < get_n(); ++i) {
                if (std::abs(eigvecs(i, k)) > vmax) {
                    vmax = std::abs(eigvecs(i, k));
                    imax = i;
                }
            }
            if (eigvecs(imax, k) < 0.0)
                eigvecs.col(k) = -eigvecs.col(k);
        }
    }

    return solver.eigenvalues();  // [a.u.], posortowane
}