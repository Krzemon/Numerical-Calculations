#pragma once
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include "functions.hpp"

namespace fs = std::filesystem;

/** * @brief Przygotowuje plik do zapisu danych. * * 
    * @param filename Nazwa pliku wewnatrz katalogu data * 
    * @param outFile Referencja do strumienia, ktory zostanie otwarty 
    */ 
inline void prepareDataFile(fs::path& dataFile, std::ofstream& outFile, const std::string& filename = "out.txt") { 
    fs::path projectRoot = PROJECT_ROOT_DIR; // makro przekazane z CMake 
    fs::path dataDir = projectRoot / "data"; 
    fs::path plotsDir = projectRoot / "plots"; 
    fs::create_directories(dataDir); 
    fs::create_directories(plotsDir); 
    dataFile = dataDir / filename; 
    outFile.open(dataFile);
    if (!outFile) { 
        std::cerr << "Nie mozna otworzyc pliku: " << dataFile << "\n"; 
        std::exit(1); 
    } 
}

inline void write_u_of_x_file(const fs::path& path, int M,
                              const std::vector<double>& xnodes,
                              const std::vector<std::vector<double>>& evecs,
                              int num_modes=5) {
    std::ofstream f(path);
    if(!f){ std::cerr<<"Cannot open "<<path<<"\n"; return; }

    f<<"# x";
    int modes = std::min(num_modes,(int)evecs.size());
    for(int mu=0; mu<modes; ++mu) f<<" u"<<mu;
    f<<"\n";

    double xa = xnodes.front();
    double xb = xnodes.back();
    double dx = 0.001; // duża liczba punktów, żeby było gładko

    for(double x=xa; x<=xb+1e-12; x+=dx){
        f << std::setprecision(12) << x;

        // znajdź element, w którym leży x
        int m_found = -1;
        for(int m=0;m<M;m++){
            double xa_m = xnodes[2*m];
            double xb_m = xnodes[2*m+2];
            if((x>=xa_m && x<xb_m)||(m==M-1 && x<=xb_m+1e-12)){
                m_found = m;
                break;
            }
        }

        if(m_found < 0){
            for(int mu=0; mu<modes; mu++) f<<" 0";
            f<<"\n";
            continue;
        }

        double xa_m = xnodes[2*m_found];
        double xb_m = xnodes[2*m_found+2];
        double xi = xi_of_x(xa_m, xb_m, x);

        for(int mu=0; mu<modes; mu++){
            double u = 0.0;
            for(int i=0;i<3;i++){
                int k = global_index(m_found,i);
                u += evecs[mu][k] * phi(i, xi);
            }
            f<<" "<<std::setprecision(12)<<u;
        }
        f<<"\n";
    }
    f.close();
}