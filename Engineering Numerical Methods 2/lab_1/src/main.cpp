#include <cmath>
#include <iostream>
#include "config.hpp"
#include "solution.cpp"

int main() {
        
// ------------------ Zadanie 1 ------------------
    fs::path zad_1;
    std::ofstream s1;
    prepareDataFile(zad_1, s1, "zad_1.dat");
    calc_zad_1(s1); // l = 0
    std::cout << "Zapisano wyniki do: " << zad_1 << "\n";

// ------------------ Zadanie 2 ------------------
    fs::path zad_2;
    std::ofstream s2;
    prepareDataFile(zad_2, s2, "zad_2.dat");
    calc_zad_2(s2); // l = 0
    std::cout << "Zapisano wyniki do: " << zad_2 << "\n";
    
// ------------------ Zadanie 3 ------------------
    fs::path zad_3_0;
    std::ofstream s30;
    prepareDataFile(zad_3_0, s30, "zad_3_0.dat");
    calc_zad_3(s30, 0); // l = 0
    std::cout << "Zapisano wyniki do: " << zad_3_0 << "\n";

    fs::path zad_3_1;
    std::ofstream s31;
    prepareDataFile(zad_3_1, s31, "zad_3_1.dat");
    calc_zad_3(s31, 1); // l = 1
    std::cout << "Zapisano wyniki do: " << zad_3_1 << "\n";

// ------------------ Zadanie 4 ------------------
    fs::path zad_4_0;
    std::ofstream s40;
    prepareDataFile(zad_4_0, s40, "zad_4_0.dat");
    calc_zad_4_1(s40, 1); // l = 1 
    std::cout << "Zapisano wyniki do: " << zad_4_0 << "\n";
    
    calc_zad_4_2(1); // l = 1

    fs::path zad_4_2;
    std::ofstream s42;
    prepareDataFile(zad_4_2, s42, "zad_4_2.dat");
    calc_zad_4_3(s42, 1); // l = 1 
    std::cout << "Zapisano wyniki do: " << zad_4_2 << "\n";

    return 0;
}



