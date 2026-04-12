#include "config.hpp"
#include "QuantumDot.hpp"
#include <iostream>
 
int main()
{
    // ----------
    // ZADANIE 1 
    // ----------
    {
        std::cout << "Zad 1: Funkcje bazowe dla siatki 9x9, dx = 2 nm, mapy dla k = 0, 8, 9\n";
        QuantumDot qdot(9, 2.0);
        for (int k : {0, 8, 9}) {
            auto f = prepareDataFile("basis_k" + std::to_string(k) + ".dat");
            qdot.write_basis(f.out, k);
            std::cout << "  " << f.path << '\n';
        }
    }
 
    // ------------
    // ZADANIE 2+3
    // ------------
    {
        std::cout << "\nZad 2: Złożenie macierzy H i S i rozwiązanie uogólnionego problemu własnego\n";

        QuantumDot qdot(9, 1.0);
        qdot.solve();
 
        auto fH = prepareDataFile("matrix_H.dat");
        auto fS = prepareDataFile("matrix_S.dat");
        qdot.write_matrix(fH.out, qdot.H, "H");
        qdot.write_matrix(fS.out, qdot.S, "S");
        std::cout << "  " << fH.path << '\n';
        std::cout << "  " << fS.path << '\n';
        
        auto fe = prepareDataFile("energies.dat");
        qdot.write_energies(fe.out, 10);
        std::cout << "  " << fe.path << '\n';

        std::cout << "\nZad 3: Mapy |psi|^2 dla 6 najnizszych, dx = 1 nm\n";

        for (int s = 0; s < 6; ++s) {
            auto f = prepareDataFile("WF_state" + std::to_string(s) + ".dat");
            qdot.write_wavefunction(f.out, s);
            std::cout << "  " << f.path << '\n';
        }
    }
 
    // ----------
    // ZADANIE 4
    // ----------
    {
        std::cout << "\nZad 4: Energie 10 stanow w funkcji hbar*wx in [0, 500] meV, wy = 200 meV\n";
        constexpr int NS = 10;
        constexpr int N_steps = 501;
        constexpr double wx1 = 10.0;
        constexpr double wx2 = 500.0;
        constexpr double dwx = (wx2 - wx1) / (N_steps - 1);
 
        auto f = prepareDataFile("energies_omega.dat");
        QuantumDot::write_energy_vs_omega_header(f.out, NS);
 
        for (int i = 0; i < N_steps; ++i) 
        {
            double hbwx_meV = wx1 + i * dwx;
            double wx = (hbwx_meV * 1e-3) / E_h;
            QuantumDot qdot(9, 1.0, m_effect, wx, omega_y);
            qdot.solve();
            qdot.write_energy_vs_omega_row(f.out, hbwx_meV, NS);
        }
        std::cout << "  " << f.path << '\n';
    }
 
    // ----------
    // ZADANIE 5
    // ----------
    {
        std::cout << "\nZad 5: Dobór wy = 400 meV tak by 5 najnizszych stanów bylo wzbudzonych tylko w x (ny=0)\n";
        QuantumDot qdot(9, 1.0, m_effect, omega_x, omega_y_zad5);
        qdot.solve();
 
        auto fe = prepareDataFile("energies_zad5.dat");
        qdot.write_energies(fe.out, 10);
        std::cout << "  " << fe.path << '\n';
 
        for (int s = 0; s < 6; ++s) 
        {
            auto f = prepareDataFile("WF_state" + std::to_string(s) + "_wy.dat");
            qdot.write_wavefunction(f.out, s);
            std::cout << "  " << f.path << '\n';
        }
    }
 
    std::cout << std::endl;
    return 0;
}