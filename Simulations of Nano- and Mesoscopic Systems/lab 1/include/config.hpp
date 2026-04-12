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
    fs::path dataDir  = projectRoot / "data";
    fs::path plotsDir = projectRoot / "plots";
    fs::create_directories(dataDir);
    fs::create_directories(plotsDir);

    fs::path zad_1_path = plotsDir / "zad_1";
    fs::path zad_2_path = plotsDir / "zad_2";
    fs::path zad_3_path = plotsDir / "zad_3";
    fs::path zad_4_path = plotsDir / "zad_4";
    fs::path zad_5_path = plotsDir / "zad_5";
    fs::create_directories(zad_1_path);
    fs::create_directories(zad_2_path);
    fs::create_directories(zad_3_path);
    fs::create_directories(zad_4_path);
    fs::create_directories(zad_5_path);
 
    FileHandle file;
    file.path = dataDir / filename;
    file.out.open(file.path);
 
    if (!file.out.is_open())
        throw std::runtime_error("Nie można otworzyć pliku: " + file.path.string());
 
    return file;
}
 