#include "solution.hpp"
#include "config.hpp"

int main() {
std::cout << "------------------ZADANIE 1------------------\n";
    auto file_1 = prepareDataFile("zad_1.dat");
    calc_zad_1(file_1.out);
    std::cout << "Zapisano: " << file_1.path << "\n";

std::cout << "------------------ZADANIE 2------------------\n";
    auto file_2 = prepareDataFile("zad_2.dat");
    calc_zad_2(file_2.out);
    std::cout << "Zapisano: " << file_2.path << "\n";

std::cout << "------------------ZADANIE 3a------------------\n";
    auto file_3a = prepareDataFile("zad_3_matrix.dat");
    calc_zad_3a(file_3a.out);
    std::cout << "Zapisano: " << file_3a.path << "\n";

std::cout << "------------------ZADANIE 3b------------------\n";
    auto file_3b = prepareDataFile("zad_3_vector.dat");
    calc_zad_3b(file_3b.out);
    std::cout << "Zapisano: " << file_3b.path << "\n";

std::cout << "------------------ZADANIE 4------------------\n";
    auto file_4 = prepareDataFile("zad_4.dat");
    calc_zad_4(file_4.out);
    std::cout << "Zapisano: " << file_4.path << "\n";

// std::cout << "------------------ZADANIE 5------------------\n";
//     auto file_5 = prepareDataFile("zad_5.dat");
//     calc_zad_5(file_5.out);
//     std::cout << "Zapisano: " << file_5.path << "\n";
    

    std::cout << std::endl;
    return 0;
}