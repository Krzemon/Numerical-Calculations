#ifndef POTENTIALSOLVER_H
#define POTENTIALSOLVER_H

#include "Grid.h"

class PotentialSolver {
public:
    PotentialSolver(Grid& grid, double tolerance, int maxIterations, double omega);

    void solve();

private:
    Grid& grid;
    double tolerance;
    int maxIterations;
    double omega;

    double calculateResidual(int x, int y) const;
};

#endif // POTENTIALSOLVER_H
