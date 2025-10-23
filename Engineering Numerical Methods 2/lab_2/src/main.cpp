#include "config.hpp"
#include "solution.cpp"

int main() {
        
std::cout << "------------------ZADANIE 1------------------\n";
    check_bc_for_base_type(1);

std::cout << "------------------ZADANIE 2------------------\n";

    fs::path zad_2;
    std::ofstream s2;
    prepareDataFile(zad_2, s2, "zad_2.dat");
    calc_collocation(s2, 1);
    std::cout << "Zapisano wyniki do: " << zad_2 << "\n";

std::cout << "------------------ZADANIE 3------------------\n";

    fs::path zad_3;
    std::ofstream s3;
    prepareDataFile(zad_3, s3, "zad_3.dat");
    calc_least_squares(s3, 1);
    std::cout << "Zapisano wyniki do: " << zad_3 << "\n";

std::cout << "------------------ZADANIE 4------------------\n";

    fs::path zad_4;
    std::ofstream s4;
    prepareDataFile(zad_4, s4, "zad_4.dat");
    calc_galerkin(s4, 1);
    std::cout << "Zapisano wyniki do: " << zad_4 << "\n";

std::cout << "------------------ZADANIE 5------------------\n";

    fs::path zad_5_1;
    std::ofstream s5_1;
    prepareDataFile(zad_5_1, s5_1, "zad_5_1.dat");
    calc_collocation(s5_1, 2);
    std::cout << "Zapisano wyniki do: " << zad_5_1 << "\n\n";

    fs::path zad_5_2;
    std::ofstream s5_2;
    prepareDataFile(zad_5_2, s5_2, "zad_5_2.dat");
    calc_least_squares(s5_2, 2);
    std::cout << "Zapisano wyniki do: " << zad_5_2 << "\n\n";

    fs::path zad_5_3;
    std::ofstream s5_3;
    prepareDataFile(zad_5_3, s5_3, "zad_5_3.dat");
    calc_galerkin(s5_3, 2);
    std::cout << "Zapisano wyniki do: " << zad_5_3 << "\n\n";

return 0;
}



