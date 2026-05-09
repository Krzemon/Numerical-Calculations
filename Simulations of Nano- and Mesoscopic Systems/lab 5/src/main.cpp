#include "config.hpp"
#include "params.hpp"
#include "QPC.hpp"

int main()
{
    stream_config(std::cout, 6, false);

    /**
     * @brief Zad 1: relacja dyspersji E(kx) dla jednorodnego kanału
     */
    ex_counter();
    {
        Grid g {
            -nm_to_au(QPCParams::y_half_nm),  nm_to_au(QPCParams::y_half_nm),
            -nm_to_au(QPCParams::y_half_nm),  nm_to_au(QPCParams::y_half_nm),
            QPCParams::Nx_small, QPCParams::Ny_small, QPCParams::m_eff
        };
        QPCSolver qpc(g);
        qpc.dispersion("zad1_dispersion.dat", 500);
    }

    /**
     * @brief Zad 2: mody propagujące się dla energii E [au]
     */
    ex_counter();
    {
        Grid g {
            -nm_to_au(QPCParams::y_half_nm),  nm_to_au(QPCParams::y_half_nm),
            -nm_to_au(QPCParams::y_half_nm),  nm_to_au(QPCParams::y_half_nm),
            QPCParams::Nx_small, QPCParams::Ny_small, QPCParams::m_eff
        };
        QPCSolver qpc(g);
        qpc.modes(eV_to_au(QPCParams::E1_eV), "zad2_modes_E1", "zad2_u_E1_");
        qpc.modes(eV_to_au(QPCParams::E2_eV), "zad2_modes_E2", "zad2_u_E2_");
    }

    /**
     * @brief Zad 3: test T+R = liczba modów
     */
    ex_counter();
    {
        Grid g {
            -nm_to_au(QPCParams::y_half_nm),  nm_to_au(QPCParams::y_half_nm),
            -nm_to_au(QPCParams::y_half_nm),  nm_to_au(QPCParams::y_half_nm),
            QPCParams::Nx_test, QPCParams::Ny_small, QPCParams::m_eff
        };
        QPCSolver qpc(g);
        
        auto fh = prepareDataFile("zad3_test_TR.dat");
        stream_config(fh.out);
        fh.out << "# E[eV]  T  R  T+R\n";
        
        for (double E_eV : { QPCParams::E1_eV, QPCParams::E2_eV }) {
            auto [T, R] = qpc.transmission(eV_to_au(E_eV), 0.0, "zad3_wf_" + std::to_string(static_cast<int>(E_eV * 1000)));
            std::cout << "E=" << E_eV << "eV  T=" << T << "  R=" << R << "  T+R=" << T+R << "\n";
            fh.out << E_eV << " " << T << " " << R << " " << T+R << "\n";
        }
        
        std::cout << " [saved]: " << fh.path << "\n";
    }

    /**
     * @brief Zad 4: mapa V i przewodność G(Vgates)
     */
    ex_counter();
    {
        Grid g {
            -nm_to_au(QPCParams::x_half_nm), nm_to_au(QPCParams::x_half_nm),
            -nm_to_au(QPCParams::y_half4_nm), nm_to_au(QPCParams::y_half4_nm),
             QPCParams::Nx_large, QPCParams::Ny_large, QPCParams::m_eff
        };
        QPCSolver qpc(g);
        qpc.conductance_vs_gate(eV_to_au(QPCParams::E4_eV),
                                QPCParams::Vg_min_eV, QPCParams::Vg_max_eV,
                                QPCParams::Vg_steps,
                                "zad4_conductance", "zad4_potential");
    }

    /**
     * @brief Zad 5: mapy gęstości |psi|^2 dla G≈0, G0, 2G0
     */
    ex_counter();
    {
        Grid g {
            -nm_to_au(QPCParams::x_half_nm), nm_to_au(QPCParams::x_half_nm),
            -nm_to_au(QPCParams::y_half4_nm), nm_to_au(QPCParams::y_half4_nm),
             QPCParams::Nx_large, QPCParams::Ny_large, QPCParams::m_eff
        };
        QPCSolver qpc(g);
        
        std::vector<std::pair<double, std::string>> cases = {
            {-1.3, "-130"},
            {-1.1, "-110"},
            {-0.85, "-85"}
        };
        
        for (const auto& [Vg_eV, label] : cases) {
            auto [T, R] = qpc.transmission(eV_to_au(QPCParams::E4_eV), eV_to_au(Vg_eV), "zad5_wf_" + label);
            std::cout << "Vg=" << Vg_eV << "eV  T=" << T << "  R=" << R << "  T+R=" << T+R << "\n";
        }
    }

    std::cout << "-----------------------------------------\n";

    return 0;
}