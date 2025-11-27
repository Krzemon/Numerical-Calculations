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
    auto file_3 = prepareDataFile("zad_3.dat");
    calc_zad_3(file_3.out);
    std::cout << "Zapisano: " << file_3.path << "\n";

std::cout << "------------------ZADANIE 4------------------\n";
    auto file_4 = prepareDataFile("zad_4.dat");
    calc_zad_4(file_4.out);
    std::cout << "Zapisano: " << file_4.path << "\n";

std::cout << "------------------ZADANIE 5------------------\n";
    calc_zad_5();

std::cout << "------------------ZADANIE 6------------------\n";
    calc_zad_6();

std::cout << "------------------ZADANIE 7------------------\n";
    auto file_7 = prepareDataFile("zad_7.dat");
    calc_zad_7(file_7.out);
    std::cout << "Zapisano: " << file_7.path << "\n";

std::cout << "------------------ZADANIE 8------------------\n";
    auto file_8_3 = prepareDataFile("zad_8_3.dat");
    calc_zad_8(file_8_3.out, 3);
    std::cout << "Zapisano: " << file_8_3.path << "\n";

    auto file_8_10 = prepareDataFile("zad_8_10.dat");
    calc_zad_8(file_8_10.out, 10);
    std::cout << "Zapisano: " << file_8_10.path << "\n";

    std::cout << std::endl;
    return 0;
}