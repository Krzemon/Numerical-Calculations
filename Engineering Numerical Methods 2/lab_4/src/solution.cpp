#include "solution.hpp"

void calc_zad_1() {
    std::cout << "Parametry zdefiniowane w pliku src/solution.cpp:\n";
}

void calc_zad_2(std::ofstream& out) {
    make_nodes(); // cooking global_nodes & local_nodes
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

void calc_zad_3() {

}

// void calc_zad_4() {
// }

// void calc_zad_5() {
// }

// void calc_zad_6() {
// }

// void calc_zad_7() {
    // solve_Ac_b(S, c, F, N);
// }

// void calc_zad_8() {
// }