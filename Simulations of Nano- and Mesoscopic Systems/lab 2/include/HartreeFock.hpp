#pragma once

#include <vector>
#include <ostream>

#include "config.hpp"
#include "params.hpp"
 
class HartreeFock {
public:
    using Grid = std::vector<double>;

    /**
     * @brief Konstruktor klasy HartreeFock z parametrami domyślnymi
     */
    explicit HartreeFock(double a_nm_in = a_nm, int N_in = n_nodes, double tol_IT_in = tolerance, double tol_HF_in = tolerance);

    /**
     * @brief Operator wypisywania danych obiektu
     */
    friend std::ostream& operator<<(std::ostream& os, const HartreeFock& hf);
 
    /**
     * @brief Uruchamia symulację metody Hatree-Focka i oblicza energię stanu podstawowego
     */
    void run();

    /**
     * @brief Zapis energii HF (zewnętrzne iteracje) do strumienia
     */
    void write_energies(std::ostream& os) const;
 
    /**
     * @brief Zapis gęstości kwadratu funkcji falowej do strumienia
     */
    void write_density(std::ostream& os, bool is_nanometers = true) const;

    /// GETTERY --------------------------

    /**
     * @brief zwraca energię stanu podstawowego w eV
     * @addtogroup inline
     */
    inline double energy() const 
    { 
        return E_final_; 
    }
// ----------------------------------------------------------------------
private:
    int N_;
    double tol_IT_, tol_HF_;
 
    double a_au_, l_au_;
    double dx_, dt_;
    long long iter_ = 0;

    Grid Psi_;
    Grid J1_;   
    Grid V_int_;
 
    std::vector<double> E_history_;
    double E_final_ = 0.0;
  
    /**
     * @brief Inicjalizuje parametry
     */
    void initialize_parameters(double a_nm_in);
 
    /**
     * @brief Inicjalizuje Psi, J1 i V_int
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
     * @brief Normalizuje przekazaną funkcję falową do wartości 1
     */
    void normalize(Grid& psi) const;
 
    /**
     * @brief Oblicza funkcję J1
     */
    void compute_J1();

    /**
     * @brief Oblicza operator Focka w działaniu na psi |F|psi>
     */
    void apply_fock(const Grid& psi, Grid& Fpsi) const;
 
    /**
     * @brief Oblicza wartość oczekiwaną operatora F, czyli energię kinetyczno-potencjalną
     */
    double compute_epsilon(const Grid& psi, const Grid& Fpsi) const;
 
    /**
     * @brief Oblicza energię oddziaływania (potencjał średni w stanie psi)
     */
    double compute_E_int() const;
 
    /**
     * @brief Wewnętrzna pętla: imaginary time 1D
     */
    double solve_fock_imaginary_time();

     /**
     * @brief Główna pętla symulacji metody Hartree-Focka
     */
    void compute();

    /**
     * @brief Sprawdza, czy punkt (i) jest na brzegu siatki
     * @addtogroup inline
     */
    inline bool is_boundary(int i) const
    {
        return i == 0 || i == N_-1;
    }
};
