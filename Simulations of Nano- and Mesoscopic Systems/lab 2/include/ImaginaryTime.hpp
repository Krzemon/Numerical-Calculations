#pragma once

#include <vector>
#include <ostream>

#include "config.hpp"
#include "params.hpp"

class ImaginaryTime {
public:
    using Grid = std::vector<double>;

    /**
     * @brief Konstruktor klasy ImaginaryTime z parametrami domyślnymi
     */
    explicit ImaginaryTime(double a_nm_in = a_nm, int N_in = n_nodes, double tol_in = tolerance);

    /**
     * @brief Operator wypisywania danych obiektu
     */
    friend std::ostream& operator<<(std::ostream& os, const ImaginaryTime& it);

    /**
     * @brief Uruchamia symulację czasu urojonego i oblicza energię stanu podstawowego
     */
    void run();

    /**
     * @brief Zapis energii w kolejnych iteracjach czasu urojonego do pliku
     */
    void write_energies(std::ostream& os) const;

    /**
     * @brief Zapis gęstości |Psi(x1, x2)|^2 do pliku
     */
    void write_density(std::ostream& os, bool is_nanometers = true) const;

    /**
     * @brief Zwraca energię stanu podstawowego w eV
     * @addtogroup inline
     */
    inline double energy() const
    {
        return E_final_;
    }
// ----------------------------------------------------------------------
private:

    int N_;;
    double tol_;

    double a_au_, l_au_;
    double dx_, dt_;
    long long iter_ = 0;

    Grid Psi_;
    Grid V_int_;
    std::vector<double> E_history_;

    double E_final_ = 0.0; // [a.u.] energia stanu podstawowego

    /**
     * @brief Inicjalizuje parametry
     */
    void initialize_parameters(double a_nm_in);

    /**
     * @brief Inicjalizuje siatkę Psi i V_int
     */
    void initialize_grid();

    /**
     * @brief Oblicza macierz oddziaływań Coulomba V_int[i-j]
     */
    void coulomb_interaction();

    /**
     * @brief Generuje losowy stan początkowy z warunkami brzegowymi Psi=0 na krawędziach
     */
    void initialize_random();

    /**
     * @brief Oblicza H*psi dla całej siatki
     */
    void apply_hamiltonian(const Grid& psi, Grid& Hpsi) const;

    /**
     * @brief Normalizuje funkcję falową do wartości 1
     */
    void normalize(Grid& psi) const;

    /**
     * @brief Oblicza energię 
     */
    double compute_energy(const Grid& psi, const Grid& Hpsi) const;

    /**
     * @brief Główna pętla symulacji czasu urojonego
     */
    void compute();

    /**
     * @brief Zwraca indeks w jednowymiarowej tablicy dla punktu (i, j) na siatce 2D
     * @addtogroup inline
     */
    inline int index(int i, int j) const
    {
        return i * N_ + j;
    }

    /**
     * @brief Sprawdza, czy punkt (i, j) jest na brzegu siatki
     * @addtogroup inline
     */
    inline bool is_boundary(int i, int j) const
    {
        return i == 0 || i == N_-1 || j == 0 || j == N_-1;
    }

};