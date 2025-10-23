#pragma once

#include <filesystem>
#include <fstream>
#include <iostream>
#include <cstdlib>

namespace fs = std::filesystem;

/**
 * @brief Przygotowuje plik do zapisu danych.
 * 
 * @param filename Nazwa pliku wewnatrz katalogu data
 * @param outFile Referencja do strumienia, ktory zostanie otwarty
 */
inline void prepareDataFile(fs::path& dataFile, std::ofstream& outFile, const std::string& filename = "out.txt") {
    fs::path projectRoot = PROJECT_ROOT_DIR;  // makro przekazane z CMake
    fs::path dataDir = projectRoot / "data";
    fs::create_directories(dataDir);

    dataFile = dataDir / filename;
    outFile.open(dataFile);
    if (!outFile) {
        std::cerr << "Nie mozna otworzyc pliku: " << dataFile << "\n";
        std::exit(1);
    }
}
