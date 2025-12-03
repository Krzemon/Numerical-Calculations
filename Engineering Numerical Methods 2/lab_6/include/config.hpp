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

/**
 * @brief Przygotowuje plik do zapisu danych w katalogu data projektu
 */
inline FileHandle prepareDataFile(const std::string& filename = "out.txt") { 
    fs::path projectRoot = PROJECT_ROOT_DIR; // makro przekazane z CMake 
    fs::path dataDir = projectRoot / "data"; 
    fs::path plotsDir = projectRoot / "plots"; 
    fs::create_directories(dataDir); 
    fs::create_directories(plotsDir); 

    FileHandle file;
    file.path = dataDir / filename;
    file.out.open(file.path);

    if (!file.out.is_open())
        throw std::runtime_error("Nie można otworzyć pliku: " + file.path.string());
        
    return file;
}