#pragma once
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>
 
namespace fs = std::filesystem;
 
struct FileHandle {
    fs::path path;
    std::ofstream out;
};

inline bool SCIENTIFIC_NOTATION = false;

/**
 * @brief konfiguruje strumień wyjściowy do zapisu danych z określoną precyzją i formatem
 */
inline void stream_config(int precision = 9, bool scientific = SCIENTIFIC_NOTATION) {
    if (scientific)
        std::cout << std::scientific;
    else
        std::cout << std::fixed;
    std::cout << std::setprecision(precision) << std::flush;
}

/**
 * @brief Przygotowuje plik do zapisu danych w katalogu data projektu
 */
inline FileHandle prepareDataFile(const std::string& filename = "out.txt") {
    fs::path projectRoot = PROJECT_ROOT_DIR; // makro przekazane z CMake
    fs::path dataDir  = projectRoot / "data";
    fs::path plotsDir = projectRoot / "plots";
    fs::create_directories(dataDir);
    fs::create_directories(plotsDir);

    fs::path zad_1_path = plotsDir / "zad_1";
    fs::path zad_2_path = plotsDir / "zad_2_4";
    fs::path zad_3_path = plotsDir / "zad_3";
    fs::create_directories(zad_1_path);
    fs::create_directories(zad_2_path);
    fs::create_directories(zad_3_path);
 
    FileHandle file;
    file.path = dataDir / filename;
    file.out.open(file.path);
 
    if (!file.out.is_open())
        throw std::runtime_error("Nie można otworzyć pliku: " + file.path.string());
 
    return file;
}

/**
 * @brief Funkcja pomocnicza do numerowania zadań w konsoli
 */
inline void ex_counter() {
    static int num = 0;
    std::cout << "--------------- \n";
    std::cout << "-- Zadanie " << ++num << " -- \n";
    std::cout << "----------------------------------------- \n";
}

