#include <iostream>
#include "zadania.hpp"
#include "prepare.hpp"

int main() {

    std::cout << "------------------ZADANIE 1------------------\n";
    auto file_1 = prepareDataFile("zad_1.dat");
    calc_zad_1(file_1.out);
    std::cout << "Zapisano: " << file_1.path << "\n\n";

    std::cout << "------------------ZADANIE 2------------------\n";
    auto file_2 = prepareDataFile("zad_2.dat");
    calc_zad_2(file_2.out);
    std::cout << "Zapisano: " << file_2.path << "\n\n";

    std::cout << "------------------ZADANIE 3------------------\n";
    auto file_3 = prepareDataFile("zad_3.dat");
    calc_zad_3(file_3.out);
    std::cout << "Zapisano: " << file_3.path << "\n\n";

    std::cout << "------------------ZADANIE 4------------------\n";
    auto file_4 = prepareDataFile("zad_4.dat");
    calc_zad_4(file_4.out);
    std::cout << "Zapisano: " << file_4.path << "\n\n";

    std::cout << "------------------ZADANIE 5------------------\n";
    auto file_5 = prepareDataFile("zad_5.dat");
    calc_zad_5(file_5.out);
    std::cout << "Zapisano: " << file_5.path << "\n\n";

    return 0;
}