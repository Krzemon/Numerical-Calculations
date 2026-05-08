#include "TransferMatrix.hpp"
#include <cmath>

using namespace std::complex_literals;

cdouble TransferMatrix::wavevec(double E_au, double U_au, double m_au)
{
    double diff = E_au - U_au;
    if (diff > 0.0) {
        return std::sqrt(2.0 * m_au * diff) / h_bar;
    } else {
        return cdouble(0.0, std::sqrt(2.0 * m_au * std::abs(diff)) / h_bar);
    }
}

Mat2 TransferMatrix::monodromy(cdouble kn, cdouble kn1, double mn, double mn1, double zn)
{
    const cdouble ratio = (kn1 * mn) / (kn * mn1);
    const cdouble coeff_plus  = 0.5 * (1.0 + ratio);
    const cdouble coeff_minus = 0.5 * (1.0 - ratio);
    
    const cdouble exp_plus  = std::exp(1.0i * (kn1 + kn) * zn);
    const cdouble exp_minus = std::exp(1.0i * (kn1 - kn) * zn);
    
    Mat2 M{};
    M[0][0] = coeff_plus  * exp_minus;
    M[0][1] = coeff_minus / exp_plus;
    M[1][0] = coeff_minus * exp_plus;
    M[1][1] = coeff_plus  / exp_minus;
    return M;
}

Mat2 TransferMatrix::mat_mul(const Mat2& A, const Mat2& B)
{
    Mat2 C{};
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            for (int k = 0; k < 2; ++k)
                C[i][j] += A[i][k] * B[k][j];
    return C;
}

Mat2 TransferMatrix::mat_identity()
{
    Mat2 I{};
    I[0][0] = 1.0;
    I[1][1] = 1.0;
    return I;
}

std::pair<double, double> TransferMatrix::TR_from_product( const Mat2& M, cdouble k1, cdouble kN, double m1, double mN)
{
    const cdouble M11 = M[0][0];
    const cdouble M21 = M[1][0];

    const double T = std::real(kN * m1 / k1 / mN) / std::norm(M11);
    const double R = std::norm(M21) / std::norm(M11);
    return {T, R};
}
