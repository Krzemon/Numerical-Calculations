#pragma once

#include "params.hpp"
#include "config.hpp"

#include <complex>
#include <array>
#include <vector>
#include <string>

using cdouble = std::complex<double>;
using Mat2 = std::array<std::array<cdouble, 2>, 2>;

/**
 * @brief Bazowa klasa metody macierzy transferu (1D).
 */
class TransferMatrix {
public:
    TransferMatrix() = default;
    virtual ~TransferMatrix() = default;

    /**
     * @brief Oblicza liczbę falową k (możliwe wartości zespolone przy E < U).
     */
    static cdouble wavevec(double E_au, double U_au, double m_au);

    /**
     * @brief Oblicza macierz monodromii M_n łączącą obszary n i n+1.
     */
    static Mat2 monodromy(cdouble kn, cdouble kn1, double mn, double mn1, double zn);

    /**
     * @brief Mnoży dwie macierze 2x2.
     */
    static Mat2 mat_mul(const Mat2& A, const Mat2& B);

    /**
     * @brief Zwraca macierz jednostkową 2x2.
     */
    static Mat2 mat_identity();

protected:
    /**
     * @brief Oblicza T i R z produktu macierzy M_{1->N}.
     */
    static std::pair<double, double> TR_from_product(const Mat2& M, cdouble k1, cdouble kN, double m1, double mN);
};
