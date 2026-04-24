#pragma once
  
/// constants
constexpr static double h_bar = 1.0;               // [a.u.] stała Plancka zredukowana
constexpr static double E_h = 27.211;              // [eV]   energia Hartree 1[Ha]
constexpr static double bohr_ratio_nm = 0.0529;    // [nm]   promien Bohra
 
/// parametry symulacji
constexpr static double tolerance = 1e-8;          // tolerancja ~e-9

struct DotParams {
    double a_nm = 25.0;       // [nm] promien kropki kwantowej
    double d1_nm = 12.0;      // [nm] zewnętrzna krawędź bariery
    double d2_nm = 4.0;       // [nm] wewnętrzna krawędź bariery
    double V1_eV = 0.25;      // [eV] głębokość bariery zewnętrznej
    double V2_eV = 0.20;      // [eV] głębokość bariery wewnętrznej
    double m_effect = 0.067;  // [m_e] efektywna masa (GaAs)
    int n = 101;              // liczba węzłów wewnętrznych (bez brzegów)
    int N = n + 1;            // rozmiar siatki z brzegami (N = n + 1)
};

struct EvolutionParams {
    double F_kVcm = 0.08;     // [kV/cm] amplituda pola
    double omega_meV = 0.00;  // [meV] częstość hbar*wx
    double dt = 1.0;          // [a.u.] krok czasowy
    long steps = 3'000'001;   // liczba kroków czasowych
    int save_every = 10'000;  // okresowy zapis
    int cn_iter = 10;         // liczba iteracji Crank-Nicolson
};

inline double au_to_eV(double e) { return e * E_h; }
inline double eV_to_au(double e) { return e / E_h; }
inline double au_to_meV(double e) { return e * E_h * 1000.0; }
inline double meV_to_au(double e) { return e / (E_h * 1000.0); }
inline double au_to_nm(double r) { return r * bohr_ratio_nm; }
inline double nm_to_au(double r) { return r / bohr_ratio_nm; }
inline double F_kVcm_to_au(double f) { return f / E_h * bohr_ratio_nm * 1e-4; }
inline double F_au_to_kVcm(double f) { return f * E_h / bohr_ratio_nm * 1e4; }
inline double au_to_ns(double t) { return t * 2.4188843265857e-8; }
