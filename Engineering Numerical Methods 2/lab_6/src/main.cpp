#include "solution.hpp"
#include "config.hpp"
#include <iostream>

int main() {

    std::cout << "------------------ZADANIE 1------------------\n";
    auto file_1 = prepareDataFile("zad_1.dat");
    calc_zad_1(file_1.out);
    std::cout << "Zapisano: " << file_1.path << "\n";

    std::cout << "------------------ZADANIE 2------------------\n";
    auto file_2_E = prepareDataFile("zad_2_global_E.dat");
    auto file_2_O = prepareDataFile("zad_2_global_O.dat");
    auto file_2_evals = prepareDataFile("zad_2_evals.dat");
    auto file_2_modes = prepareDataFile("zad_2_modes.dat");
    calc_zad_2(file_2_E.out, file_2_O.out, file_2_evals.out, file_2_modes.out);
    std::cout << "Zapisano: " << file_2_E.path << "\n";
    std::cout << "Zapisano: " << file_2_O.path << "\n";
    std::cout << "Zapisano: " << file_2_evals.path << "\n";

    // std::cout << "------------------ZADANIE 3------------------\n";
    // auto file_3 = prepareDataFile("zad_3.dat");
    // calc_zad_3(file_3.out);
    // std::cout << "Zapisano: " << file_3.path << "\n";

    // std::cout << "------------------ZADANIE 4------------------\n";
    // auto file_4 = prepareDataFile("zad_4.dat");
    // calc_zad_4(file_4.out);
    // std::cout << "Zapisano: " << file_4.path << "\n";

    // std::cout << "------------------ZADANIE 5------------------\n";
    // auto file_5 = prepareDataFile("zad_5.dat");
    // calc_zad_5(file_5.out);
    // std::cout << "Zapisano: " << file_5.path << "\n";

    // std::cout << "------------------ZADANIE 6------------------\n";
    // auto file_6 = prepareDataFile("zad_6.dat");
    // calc_zad_6(file_6.out);
    // std::cout << "Zapisano: " << file_6.path << "\n";

    std::cout << std::endl;
    return 0;
}