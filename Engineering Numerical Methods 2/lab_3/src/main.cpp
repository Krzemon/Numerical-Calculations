#include <iostream>
#include <fstream>
#include <filesystem>
#include "functions.hpp"
#include "config.hpp"

extern void calc_zad_1(std::ofstream& out);
extern void calc_zad_2(std::ofstream& out, int M, const std::string& psi_prefix="zad_2_psi_M");
extern void calc_zad_3(std::ofstream& out, int M, const std::string& psi_prefix="zad_3_psi_M");

int main() {
std::cout << "------------------ZADANIE 1------------------\n";
    fs::path path1; std::ofstream f1;
    prepareDataFile(path1,f1,"zad_1_matrix.dat");
    calc_zad_1(f1);
    f1.close();
    std::cout << "Zapisano: " << path1 << "\n";

std::cout << "------------------ZADANIE 2------------------\n";
    std::vector<int> Ms2 = {5};
    for(int M: Ms2){
        fs::path pathE; std::ofstream fE;
        prepareDataFile(pathE, fE,"zad_2_E_M" + std::to_string(M) + ".dat");
        calc_zad_2(fE, M, "zad_2_psi_M");
        fE.close();
        std::cout << "Zapisano pliki dla M=" << M << "\n";
    }

std::cout << "------------------ZADANIE 3------------------\n";
    std::vector<int> Ms3 = {10,30};
    for(int M: Ms3){
        fs::path pathE; std::ofstream fE;
        prepareDataFile(pathE, fE,"zad_3_E_M" + std::to_string(M) + ".dat");
        calc_zad_3(fE, M, "zad_3_psi_M");
        fE.close();
        std::cout << "Zapisano pliki dla M=" << M << "\n";
    }

    return 0;
}