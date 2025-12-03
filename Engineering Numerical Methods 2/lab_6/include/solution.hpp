#pragma once
#include "functions.hpp"
#include <fstream>

void save_global_E(std::ofstream& out); // zapisuje do pliku globalna macierz E sztywnosci
void save_global_O(std::ofstream& out); // zapisuje do pliku globalna macierz O calek przekrywania 

void calc_zad_1(std::ofstream &out);
void calc_zad_2(std::ofstream &E_out, std::ofstream &O_out, std::ofstream &evals_out, std::ofstream &modes_out);
void calc_zad_3(std::ofstream &out);
void calc_zad_4(std::ofstream &out);
void calc_zad_5(std::ofstream &out);
void calc_zad_6(std::ofstream &out);