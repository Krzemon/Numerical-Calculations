#include "functions.hpp"

#define N_SIZE 5
constexpr std::array<int,N_SIZE> N = {6, 7, 8, 9, 10}; // liczba funkcji bazowych

inline void check_bc_for_base_type(int base_type) {
    std::cout << "Sprawdzanie warunkow brzegowych dla funkcji bazowej 1:\n";
    int n = 10;
    std::vector<std::function<double(double)>> v;
    if(base_type == 1) {
        v = gen_base_1(n);
    }else if (base_type == 2) {
        v = gen_base_2(n);
    }else {
        std::cerr << "Nieznany typ funkcji bazowej!\n";
        return;
    }
    check_bc(v, n, base_type);
}

inline void calc_collocation(std::ostream& file, int base_type) {
    std::cout << "Metoda kolokacji z funkcjami bazowymi " << base_type << "\n";
    std::vector<std::function<double(double)>> v;
    for(int w = 0; w < N_SIZE; ++w){
        int n = N[w];
        std::vector<double> x_nodes(n);
        double dx = (2*M_PI)/(n+1); // krok wezlow (brzegowe nie liczymy)
        if(base_type == 1) {
            v = gen_base_1(n);
        }else if (base_type == 2) {
            v = gen_base_2(n);
        }else {
            std::cerr << "Nieznany typ funkcji bazowej!\n";
            return;
        }

        for(int k = 0; k < n; ++k)
            x_nodes[k] = -M_PI + (k+1)*dx;

        std::vector<std::vector<double>> A(n, std::vector<double>(n, 0.0));
        std::vector<double> c(n, 0.0); // wspolczynniki rozwiazania
        std::vector<double> b(n, 0.0);

        for(int j = 0; j < n; ++j) {
            for(int i = 0; i < n; ++i) {
                A[j][i] = diff_u2(v[i], x_nodes[j]) -  diff_u1(v[i], x_nodes[j]);
            }
            b[j] = f(x_nodes[j]);
            // b[j] = 2*(1 - 3*x_nodes[j] + x_nodes[j]*x_nodes[j] - M_PI*M_PI) * std::exp(-x_nodes[j]);
        }

        solve_Ax_b(A, c, b, n);
        double max_res = check_max_row_residuum(A, c, b, n);
        std::cout << "Maksymalne residuum dla N = " << n << " : " << max_res << std::endl;

        // Zapis do wierszy pliku: N, x, u_d, u_n
        file << "# N = " << n << "\n";
        std::vector<double> x_vec = gen_x(-M_PI, M_PI);
        for (auto x : x_vec)
            file << x << " ";
        file << "\n";
        for (auto x : x_vec) {
            double u_d = analytical_solution(x);
            file << u_d << " ";
        }
        file << "\n";
        for (auto x : x_vec) {
            double u_n = 0.0;
            for (int i = 0; i < n; ++i)
                u_n += c[i] * v[i](x);
            file << u_n << " ";
        }
        file << "\n\n";
    }
}

inline void calc_least_squares(std::ostream& file, int base_type) {
    std::cout << "Metoda najmniejszych kwadratow z funkcjami bazowymi " << base_type << "\n";
    std::vector<std::function<double(double)>> v;
    for(int w = 0; w < N_SIZE; ++w){
        int n = N[w];
        if(base_type == 1) {
            v = gen_base_1(n);
        }else if (base_type == 2) {
            v = gen_base_2(n);
        }else {
            std::cerr << "Nieznany typ funkcji bazowej!\n";
            return;
        }

        auto Lvi = [&](int k, double x, double dx = delta_x) -> double {
            return diff_u2(v[k], x, dx) - diff_u1(v[k], x, dx);
        };

        auto integrate_Lvk_Lvi = [=](int k, int i, double dx = delta_x) -> double {
            size_t n_int = 40;
            double integral, xi, wi;
            gsl_integration_glfixed_table *tab= gsl_integration_glfixed_table_alloc(n_int);
            integral=0.;
            for(size_t j = 0; j < n_int; ++j){
                gsl_integration_glfixed_point(-M_PI, M_PI, j, &xi, &wi, tab);
                integral += wi * Lvi(k,xi,dx) * Lvi(i,xi,dx);
            }
            gsl_integration_glfixed_table_free(tab);
            return integral;
        };
        
        auto integrate_f_Lvk = [=](int k, double dx = delta_x) -> double {
            size_t n_int = 40;
            double integral, xi, wi;
            gsl_integration_glfixed_table *tab= gsl_integration_glfixed_table_alloc(n_int);
            integral=0.;
            for(size_t j = 0; j < n_int; ++j){
                gsl_integration_glfixed_point(-M_PI, M_PI, j, &xi, &wi, tab);
                integral += wi * f(xi) * Lvi(k,xi,dx);
            }
            gsl_integration_glfixed_table_free(tab);
            return integral;
        };

        std::vector<std::vector<double>> A(n, std::vector<double>(n, 0.0));
        std::vector<double> c(n, 0.0); // wspolczynniki rozwiazania
        std::vector<double> b(n, 0.0);

        for (int k = 0; k < n; ++k) {
            for (int i = 0; i < n; ++i) {
                A[k][i] = integrate_Lvk_Lvi(k, i);
            }
            b[k] = integrate_f_Lvk(k);
        } 

        solve_Ax_b(A, c, b, n);
        double max_res = check_max_row_residuum(A, c, b, n);
        std::cout << "Maksymalne residuum dla N = " << n << " : " << max_res << std::endl;

        // Zapis do wierszy pliku: N, x, u_d, u_n
        file << "# N = " << n << "\n";
        std::vector<double> x_vec = gen_x(-M_PI, M_PI);
        for (auto x : x_vec)
            file << x << " ";
        file << "\n";
        for (auto x : x_vec) {
            double u_d = analytical_solution(x);
            file << u_d << " ";
        }
        file << "\n";
        for (auto x : x_vec) {
            double u_n = 0.0;
            for (int i = 0; i < n; ++i)
                u_n += c[i] * v[i](x);
            file << u_n << " ";
        }
        file << "\n\n";
    }
}

inline void calc_galerkin(std::ostream& file, int base_type) {
    std::cout << "Metoda Galerkina z funkcjami bazowymi " << base_type << "\n";
    std::vector<std::function<double(double)>> v;
    for(int w = 0; w < N_SIZE; ++w){
        int n = N[w];
        if(base_type == 1) {
            v = gen_base_1(n);
        }else if (base_type == 2) {
            v = gen_base_2(n);
        }else {
            std::cerr << "Nieznany typ funkcji bazowej!\n";
            return;
        }
        
        auto Lvi = [&](int k, double x, double dx = delta_x) -> double {
            return diff_u2(v[k], x, dx) - diff_u1(v[k], x, dx);
        };

        auto integrate_vk_Lvi = [=](int k, int i, double dx = delta_x) -> double {
            size_t n_int = 40;
            double integral, xi, wi;
            gsl_integration_glfixed_table *tab= gsl_integration_glfixed_table_alloc(n_int);
            integral=0.;
            for(size_t j = 0; j < n_int; ++j){
                gsl_integration_glfixed_point(-M_PI, M_PI, j, &xi, &wi, tab);
                integral += wi * v[k](xi) * Lvi(i,xi,dx);
            }
            gsl_integration_glfixed_table_free(tab);
            return integral;
        };
        
        auto integrate_f_vk = [=](int k, double dx = delta_x) -> double {
            size_t n_int = 40;
            double integral, xi, wi;
            gsl_integration_glfixed_table *tab= gsl_integration_glfixed_table_alloc(n_int);
            integral=0.;
            for(size_t j = 0; j < n_int; ++j){
                gsl_integration_glfixed_point(-M_PI, M_PI, j, &xi, &wi, tab);
                integral += wi * f(xi) * v[k](xi);
            }
            gsl_integration_glfixed_table_free(tab);
            return integral;
        };

        std::vector<std::vector<double>> A(n, std::vector<double>(n, 0.0));
        std::vector<double> c(n, 0.0); // wspolczynniki rozwiazania
        std::vector<double> b(n, 0.0);

        for (int k = 0; k < n; ++k) {
            for (int i = 0; i < n; ++i) {
                A[k][i] = integrate_vk_Lvi(k, i);
            }
            b[k] = integrate_f_vk(k);
        } 

        solve_Ax_b(A, c, b, n);
        double max_res = check_max_row_residuum(A, c, b, n);
        std::cout << "Maksymalne residuum dla N = " << n << " : " << max_res << std::endl;

        // Zapis do wierszy pliku: N, x, u_d, u_n
        file << "# N = " << n << "\n";
        std::vector<double> x_vec = gen_x(-M_PI, M_PI);
        for (auto x : x_vec)
            file << x << " ";
        file << "\n";
        for (auto x : x_vec) {
            double u_d = analytical_solution(x);
            file << u_d << " ";
        }
        file << "\n";
        for (auto x : x_vec) {
            double u_n = 0.0;
            for (int i = 0; i < n; ++i)
                u_n += c[i] * v[i](x);
            file << u_n << " ";
        }
        file << "\n\n";
    }
}

