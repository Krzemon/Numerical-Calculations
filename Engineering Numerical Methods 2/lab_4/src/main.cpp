#include "solution.hpp"
#include "config.hpp"

int main() {
std::cout << "------------------ZADANIE 1------------------\n";
    calc_zad_1();

std::cout << "------------------ZADANIE 2------------------\n";
    auto file_2 = prepareDataFile("zad_2.dat");
    calc_zad_2(file_2.out);
    std::cout << "Zapisano: " << file_2.path << "\n";

std::cout << "------------------ZADANIE 3------------------\n";
    // auto file_2 = prepareDataFile("zad_2.dat");
    // calc_zad_3();
    // calc_zad_2(file_2.out);
    // std::cout << "Zapisano: " << file_2.path << "\n";

std::cout << "------------------ZADANIE 7------------------\n";
    // auto file_7 = prepareDataFile("zad_1.dat");
    // calc_zad_1(file_1.out);
    // std::cout << "Zapisano: " << file_1.path << "\n";

    assemble_global_matrices_MES_2D();
    solve_Ac_b(S, c, F, N_max);


    std::cout << std::endl;
    return 0;
}