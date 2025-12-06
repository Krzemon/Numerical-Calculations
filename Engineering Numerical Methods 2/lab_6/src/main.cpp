#include "solution.hpp"
#include "config.hpp"
#include <iostream>

int main() {

    std::cout << "------------------ZADANIE 1------------------\n";
    auto file_1 = prepareDataFile("zad_1.dat");
    calc_zad_1(file_1.out);
    std::cout << "Zapisano: " << file_1.path << "\n";

    // std::cout << "------------------ZADANIE 2------------------\n";
    // auto file_2_E = prepareDataFile("zad_2_global_E.dat");
    // auto file_2_O = prepareDataFile("zad_2_global_O.dat");
    // auto file_2_evals = prepareDataFile("zad_2_evals.dat");
    // calc_zad_2(file_2_E.out, file_2_O.out, file_2_evals.out);
    // std::cout << "Zapisano: " << file_2_E.path << "\n";
    // std::cout << "Zapisano: " << file_2_O.path << "\n";
    // std::cout << "Zapisano: " << file_2_evals.path << "\n";

    // std::cout << "------------------ZADANIE 3------------------\n";
    // calc_zad_3();
    // std::cout << "Znormalizowano wektory c_2 oraz c_3"<< "\n";

    // std::cout << "------------------ZADANIE 4------------------\n";
    // calc_zad_4();
    // std::cout << "Wygenerowano warunki poczatkowe y_0 oraz v_0"<< "\n";

    // std::cout << "------------------ZADANIE 5------------------\n";
    // auto file_5 = prepareDataFile("zad_5.dat");
    // calc_zad_5(file_5.out, true);
    // std::cout << "Zapisano: " << file_5.path << "\n";

    // std::cout << "------------------ZADANIE 6------------------\n";
    // auto file_6_1 = prepareDataFile("zad_6_1.dat");
    // auto file_6_2 = prepareDataFile("zad_6_2.dat");
    // auto file_6_3  = prepareDataFile("zad_6_3.dat");
    // auto file_6_4  = prepareDataFile("zad_6_4.dat");
    // auto file_6_5  = prepareDataFile("zad_6_5.dat");
    // calc_zad_6(file_6_1.out, file_6_2.out, file_6_3.out, file_6_4.out, file_6_5.out);
    // std::cout << "Zapisano: " << file_6_1.path << "\n";
    // std::cout << "Zapisano: " << file_6_2.path << "\n";
    // std::cout << "Zapisano: " << file_6_3.path << "\n";
    // std::cout << "Zapisano: " << file_6_4.path << "\n";
    // std::cout << "Zapisano: " << file_6_5.path << "\n";

    // std::cout << "------------------ GIF ------------------\n";
    // auto gif_file = prepareDataFile("GIF.dat");
    // write_gif_snapshots(gif_file.out);
    // std::cout << "Zapisano: " << gif_file.path << "\n";

    std::cout << std::endl;
    return 0;
}