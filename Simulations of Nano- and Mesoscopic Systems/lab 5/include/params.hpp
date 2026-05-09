#pragma once

constexpr static double h_bar = 1.0;
constexpr static double E_h = 27.211;
constexpr static double bohr_ratio_nm = 0.0529177;

struct QPCParams {
    static constexpr double m_eff = 0.017;     // masa efektywna InSb

    // Zadanie 1-2
    static constexpr double y_half_nm = 10.0;  // połowa szerokości kanału
    static constexpr int    Ny_small  = 19;    // liczba węzłów w y
    static constexpr int    Nx_small  = 3;     // liczba węzłów w x (dla zad 1-2)
    static constexpr int    Nx_test   = 9;     // liczba węzłów w x (dla zad 3)

    // Zadanie 2
    static constexpr double E1_eV = 0.2; // energia dla 1 modu (200 meV)
    static constexpr double E2_eV = 0.4; // energia dla 2 modów (400 meV)

    // Zadanie 3
    static constexpr double sigma_nm    = 300.0;  // sigma dla potencjału bramek
    static constexpr double V_prefactor = -0.035; // współczynnik przed exp w V_QPC

    // Zadanie 4-5
    static constexpr double x_half_nm  = 500.0; // połowa długości QPC
    static constexpr double y_half4_nm = 350.0; // połowa szerokości QPC
    static constexpr int    Nx_large   = 49;    // liczba węzłów w x
    static constexpr int    Ny_large   = 34;    // liczba węzłów w y
    static constexpr double E4_eV      = 0.015; // energia dla obliczeń T i R
    static constexpr double Vg_min_eV  = -1.3;  // min napięcie bramek
    static constexpr double Vg_max_eV  = -0.7;  // max napięcie bramek
    static constexpr int    Vg_steps   = 61;    // liczba kroków Vgates (krok 0.01 eV)
};

inline double au_to_eV(double e) { return e * E_h; }
inline double eV_to_au(double e) { return e / E_h; }
inline double au_to_meV(double e) { return e * E_h * 1000.0; }
inline double meV_to_au(double e) { return e / (E_h * 1000.0); }
inline double au_to_nm(double r) { return r * bohr_ratio_nm; }
inline double nm_to_au(double r) { return r / bohr_ratio_nm; }
