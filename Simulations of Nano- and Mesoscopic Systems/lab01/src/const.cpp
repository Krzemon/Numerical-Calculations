#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

constexpr double Ha_to_meV_CF = 27211.6 ; // Hartree to meV conversion factor
constexpr double ab_to_nm_CF = 0.05292; // Angstrom to nm conversion factor
E_h = 27211.6
typedef std::pair<double, double> vec2d_t;
typedef std::pair<size_t, size_t> vec2ui_t;

inline double Ha_to_meV(double E)
{
    return E * Ha_to_meV_CF;
}

inline double meV_to_Ha(double E)
{
    return E / Ha_to_meV_CF;
}

inline double ab_to_nm(double L)
{
    return L * ab_to_nm_CF;
}

inline double nm_to_ab(double L)
{
    return L / ab_to_nm_CF;
}

#endif // CONSTANTS_HPP