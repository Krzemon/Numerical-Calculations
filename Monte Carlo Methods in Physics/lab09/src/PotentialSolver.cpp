#include "PotentialSolver.h"
#include <cmath>
#include <iostream>

PotentialSolver::PotentialSolver(Grid& grid, double tolerance, int maxIterations, double omega)
    : grid(grid), tolerance(tolerance), maxIterations(maxIterations), omega(omega) {}

void PotentialSolver::solve() {
    int it = 0;
    double residual = 1.0;
    double previousResidual = 0.0;

    while (residual > tolerance && it < maxIterations) {
        residual = 0.0;

        for (int j = 1; j < grid.getSizeY() - 1; ++j) {
            for (int i = 1; i < grid.getSizeX() - 1; ++i) {
                double newPotential = (1 - omega) * grid.getPotential(i, j)
                    + omega * 0.25 * (grid.getPotential(i+1, j) + grid.getPotential(i-1, j)
                                      + grid.getPotential(i, j+1) + grid.getPotential(i, j-1));

                // Calculate residual for the current point
                double pointResidual = std::abs(newPotential - grid.getPotential(i, j));
                residual += pointResidual * pointResidual;

                // Update the potential
                grid.setPotential(i, j, newPotential);
            }
        }

        residual = sqrt(residual);
        ++it;

        // Check if the residual is increasing, which indicates divergence
        if (it > 1 && residual > previousResidual) {
            std::cerr << "SOR method is diverging. Consider reducing omega or increasing iterations/tolerance." << std::endl;
            break;
        }

        previousResidual = residual;
    }

    if (it >= maxIterations) {
        std::cerr << "SOR method did not converge within max iterations." << std::endl;
    }
}

double PotentialSolver::calculateResidual(int x, int y) const {
    return pow(grid.getPotential(x, y), 2) - 0.25 * (pow(grid.getPotential(x+1, y), 2)
        + pow(grid.getPotential(x-1, y), 2) + pow(grid.getPotential(x, y+1), 2)
        + pow(grid.getPotential(x, y-1), 2));
}
