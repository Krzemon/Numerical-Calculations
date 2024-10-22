#ifndef GRID_H
#define GRID_H

#include <vector>
#include <string>

class Grid {
private:
    int nx;
    int ny;
    double delta;
    std::vector<std::vector<double>> potential;

public:
    Grid(int nx, int ny, double delta);

    void setPotential(int x, int y, double value);
    double getPotential(int x, int y) const;
    int getSizeX() const { return nx; }
    int getSizeY() const { return ny; }

    void printGrid() const;
    void savePotentialToFile(const std::string& filename) const;
};

#endif // GRID_H
