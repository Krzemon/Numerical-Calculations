#pragma once
 
// Mapowanie indeksow: k = i*n_grid + j, i = k / n_grid, j = k % n_grid
 
constexpr double h_bar = 1.0;       // [a.u.]
constexpr double E_h   = 27.211;    // [eV]   energia Hartree
constexpr double a_0   = 0.0529;    // [nm]   promien Bohra
 
constexpr double m_effect     = 0.24;            // [a.u.] masa efektywna m*
constexpr double omega_x      = 0.080 / E_h;     // [a.u.] hbar*wx = 80  meV
constexpr double omega_y      = 0.200 / E_h;     // [a.u.] hbar*wy = 200 meV
constexpr double omega_y_zad5 = 0.400 / E_h;     // [a.u.] hbar*wy = 400 meV
 
constexpr int n = 9;   // wezly na kierunek, N = n*n lacznie