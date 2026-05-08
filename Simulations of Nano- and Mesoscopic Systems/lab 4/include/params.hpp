#pragma once

// =====================================================
//  Stałe fizyczne (jednostki atomowe)
// =====================================================
constexpr static double h_bar = 1.0;
constexpr static double E_h = 27.211; // Hartree
constexpr static double bohr_ratio_nm = 0.0529177;

// =====================================================
//  Parametry geometryczne struktur
// =====================================================
struct GeometryParams {
    static constexpr double L_contact_nm = 20.0;
    static constexpr double d_barrier_nm = 5.0;
    static constexpr double d_well_nm = 3.0;
    static constexpr double U_barrier_eV = 0.27;
    static constexpr double m_GaAs = 0.063;
    static constexpr double m_AlGaAs_x = 0.3;
};

// =====================================================
//  Parametry RTD (Resonant Tunneling Diode)
// =====================================================
struct RTDParams {
    static constexpr double d_barrier_nm = 5.0;
    static constexpr double d_well_nm = 3.0;
    static constexpr double U_barrier_eV = 0.27;
    static constexpr double m_GaAs = 0.063;
    static constexpr double m_AlGaAs_x = 0.3;
    static constexpr int N_slices = 530;
    
    static constexpr double mu_eV = 0.087;
    static constexpr double T_He_K = 4.2;
    static constexpr double T_N_K = 77.0;
};

// =====================================================
//  Parametry QPC (Quantum Point Contact)
// =====================================================
struct QPCParams {
    static constexpr double m_eff = 0.063;
    static constexpr double eps = 13.6;
    static constexpr double d_nm = 3.0;
    static constexpr double W_nm = 50.0;
    static constexpr double L_nm = 100.0;
    
    static constexpr double l_frac = 0.3;
    static constexpr double r_frac = 0.7;
    static constexpr double tu_frac = 1.2;
    static constexpr double bu_frac = 0.8;
    static constexpr double tl_frac = 0.2;
    static constexpr double bl_frac = -0.2;
    
    static constexpr int N_subbands = 5;
    static constexpr int Nx = 500;
    static constexpr int Ny = 99;
    
    static constexpr double V_g_eV = 4.0;
};

// =====================================================
//  Parametry symulacji
// =====================================================
constexpr static double tolerance = 1e-9;

constexpr static auto calc_m_AlGaAs = [](double x) { return 0.063 + 0.083 * x;};

// =====================================================
//  Przeliczniki jednostek
// =====================================================
inline double au_to_eV (double e) { return e * E_h; }
inline double eV_to_au (double e) { return e / E_h; }
inline double au_to_meV(double e) { return e * E_h * 1000.0; }
inline double meV_to_au(double e) { return e / (E_h * 1000.0); }
inline double au_to_nm (double r) { return r * bohr_ratio_nm; }
inline double nm_to_au (double r) { return r / bohr_ratio_nm; }
