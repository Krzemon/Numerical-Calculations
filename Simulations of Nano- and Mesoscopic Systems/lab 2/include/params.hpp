#pragma once
  
/// constants
constexpr double h_bar = 1.0;               // [a.u.] stała Plancka zredukowana
constexpr double E_h = 27.211;              // [eV]   energia Hartree 1[Ha]
constexpr double bohr_ratio_nm = 0.0529;    // [nm]   promien Bohra
 
/// GaAs
constexpr double m_effect = 0.067;  // [a.u.] masa efektywna m*
constexpr double eps_r = 12.5;      // [a.u.] stala dielektryczna materiału

/// parametry symulacji
constexpr int n_nodes = 41;         // wezly; N = n*n lacznie
constexpr double tolerance = 1e-9;  // tolerancja ~e-9; może być za długo lepiej ~-6
constexpr double a_nm = 30.0;       // [nm] polowa szerokosci siatki [-a, a]
constexpr double l_nm = 10.0;       // [nm] szerokosc poprzeczna drutu  
