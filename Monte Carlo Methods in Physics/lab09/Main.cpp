#include "Grid.h"
#include "PotentialSolver.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>

void saveMatrixToFile(const std::vector<std::vector<double>>& matrix, const std::string& filename) {
    std::ofstream file(filename);
    if (!file) {
        std::cerr << "Cannot open file: " << filename << std::endl;
        return;
    }

    for (size_t j = 0; j < matrix[0].size(); ++j) {
        for (size_t i = 0; i < matrix.size(); ++i) {
            file << std::setw(12) << std::setprecision(6) << matrix[i][j] << " ";
        }
        file << std::endl;
    }

    file.close();
}

int main() {
    // Define parameters
    int nx = 30;
    int ny = 30;
    double delta = 0.1;
    double xmax = delta * nx;
    double ymax = delta * ny;
    double tolerance = 1e-6;
    int maxIterations = 10000;
    double omega = 1.8;

    // Create grid
    Grid grid(nx, ny, delta);

    // Set boundary conditions
    for (int i = 0; i < nx; ++i) {
        grid.setPotential(i, 0, -1.0);     // Bottom boundary
        grid.setPotential(i, ny-1, -1.0);  // Top boundary
    }
    for (int j = 0; j < ny; ++j) {
        grid.setPotential(0, j, 1.0);      // Left boundary
        grid.setPotential(nx-1, j, 1.0);   // Right boundary
    }

    // Solve using SOR method
    PotentialSolver solver(grid, tolerance, maxIterations, omega);
    solver.solve();

    // Save SOR potential to file
    std::vector<std::vector<double>> sorPotential(nx, std::vector<double>(ny));
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            sorPotential[i][j] = grid.getPotential(i, j);
        }
    }
    saveMatrixToFile(sorPotential, "../data/potential_sor.txt");

    // Example: Monte Carlo potential, deviation, absorbed chains fraction
    // Here you would implement your Monte Carlo method and save its results accordingly
    // For demonstration purposes, I'll show a placeholder for Monte Carlo results

    // 1. Monte Carlo potential (placeholder data)
    std::vector<std::vector<double>> mcPotential(nx, std::vector<double>(ny, 0.0));
    // Fill mcPotential with your Monte Carlo results
    saveMatrixToFile(mcPotential, "../data/potential_mc.txt");

    // 2. Difference between Monte Carlo and SOR potential (placeholder data)
    std::vector<std::vector<double>> potentialDiff(nx, std::vector<double>(ny, 0.0));
    // Fill potentialDiff with differences
    saveMatrixToFile(potentialDiff, "../data/potential_diff.txt");

    // 3. Standard deviation of Monte Carlo potential (placeholder data)
    std::vector<std::vector<double>> stdDeviation(nx, std::vector<double>(ny, 0.0));
    // Fill stdDeviation with standard deviation data
    saveMatrixToFile(stdDeviation, "../data/std_deviation.txt");

    // 4. Fraction of absorbed chains (placeholder data)
    std::vector<std::vector<double>> absorbedChains(nx, std::vector<double>(ny, 0.0));
    // Fill absorbedChains with fraction data
    saveMatrixToFile(absorbedChains, "../data/absorbed_chains.txt");

    std::cout << "SOR potential and Monte Carlo results saved to files." << std::endl;

    return 0;
}
