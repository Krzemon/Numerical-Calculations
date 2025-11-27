#include "solution.hpp"

void save_nodes(std::ofstream& out) {
    out << std::fixed << std::setprecision(6);
    out << "# element global_node local_node x y" << "\n";

    for (unsigned int elem_idx = 1; elem_idx <= local_nodes.size(); ++elem_idx) {
        for (unsigned int local_idx = 1; local_idx <= local_nodes[elem_idx-1].size(); ++local_idx) {
            const Node* node = local_nodes[elem_idx-1][local_idx-1];
            out << elem_idx
                << std::setw(13) << node->idx
                << std::setw(13) << local_idx
                << std::setw(13) << node->x
                << std::setw(13) << node->y
                << "\n";
        }
    }
}

void calc_zad_1() {
    std::cout << "Parametry zdefiniowane w pliku src/functions.cpp:\n";
}

void calc_zad_2(std::ofstream& out) {
    make_nodes(); // cooking global_nodes & local_nodes
    save_nodes(out);
    std::cout << "Wygenerowano siatke (elementy, węzły lokalne i globalne)\n";
}

void calc_zad_3() {
    assemble_global_matrices_MES_2D();
    std::cout << "Złożono macierz sztywności S oraz wektor obciążeń F\n";
}

void calc_zad_4() {
    border_conditions();
    std::cout << "Uwzględniono warunki brzegowe dla macierzy sztywności S oraz wektora obciążeń F\n";
}

void calc_zad_5() {
    solve_Ac_b(S, c, F, N_max);
    std::cout << "Rozwiązano układ równań Sc=F\n";
}

void calc_zad_6() {
    double a_num = functional_integral();
    std::cout << "Obliczono całkę funkcjonalną dla nx=ny=3 a_num = "<< a_num <<"\n";
}

void calc_zad_7(std::ofstream& out) {
    std::vector n_table = {5,10,15,20};
    out << std::fixed << std::setprecision(6);
    out << "# nx ny a_num" << "\n";
    for (auto n: n_table) {
        update_params(n, n);
        make_nodes();
        assemble_global_matrices_MES_2D();
        border_conditions();
        solve_Ac_b(S, c, F, N_max);
        double a_num = functional_integral();
        std::cout << "Obliczono całkę funkcjonalną dla nx=ny="<<n<<
                     " a_num = "<< a_num <<"\n";
        out << n << " "<< n << " " << a_num << "\n";
    }
}

void calc_zad_8(std::ofstream& out, int n) {

    update_params(n, n);
    make_nodes();
    assemble_global_matrices_MES_2D();
    border_conditions();
    solve_Ac_b(S, c, F, N_max);

    int n_points_x = 401;
    int n_points_y = 401;
    double dx = (x_max - x_min) / (n_points_x - 1);
    double dy = (y_max - y_min) / (n_points_y - 1);

    out << std::fixed << std::setprecision(6);
    out << "# x y u_numerical u_analytical n=" + std::to_string(n) << std::endl;

    for (int i = 0; i < n_points_x; ++i) {
        double x = x_min + i*dx;
        for (int j = 0; j < n_points_y; ++j) {
            double y = y_min + j*dy;
            out << x << " " << y << " " << u_solution(x, y) << " " << u_analytical(x, y) << std::endl;
        }
    }
    std::cout << "Zapisano dane potrzebne do wygenerowania wykresu konturowego rozwiązania numerycznego dla nx = ny =" << n << " oraz rozwiązania dokładnego.\n";
}