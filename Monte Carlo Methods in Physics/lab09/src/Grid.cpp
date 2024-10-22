#include "Grid.h"
#include <iostream>
#include <fstream>
#include <iomanip>

Grid::Grid(int nx, int ny, double delta)
    : nx(nx), ny(ny), delta(delta)
{
    potential.resize(nx, std::vector<double>(ny, 0.0));
}

void Grid::setPotential(int x, int y, double value) {
    potential[x][y] = value;
}

double Grid::getPotential(int x, int y) const {
    return potential[x][y];
}

void Grid::printGrid() const {
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            std::cout << std::setw(8) << std::setprecision(4) << potential[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

void Grid::savePotentialToFile(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file) {
        std::cerr << "Cannot open file: " << filename << std::endl;
        return;
    }

    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            file << std::setw(8) << std::setprecision(4) << potential[i][j] << " ";
        }
        file << std::endl;
    }

    file.close();
}
