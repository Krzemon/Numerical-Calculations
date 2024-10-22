#include "LJ_Potential.h"

LJ_Potential::LJ_Potential(double sigma, double eps, double rcut) {
    if (sigma <= 0.0 || eps <= 0.0 || rcut <= 0.0) {
        throw std::invalid_argument("Parametry sigma, eps i rcut muszą być dodatnie.");
    }

    this->sigma = sigma;
    this->eps = eps;
    this->rcut = rcut;

    for (double d = 0.01; d <= rcut; d += 0.1) {
        double r6 = pow(sigma / d, 6);
        double F = 24.0 * eps * r6 * (2.0 * r6 - 1) / d;
        forceValues[d] = F;
    }
}

LJ_Potential& LJ_Potential::getInstance(double sigma, double eps, double rcut) {
    static LJ_Potential instance(sigma, eps, rcut);
    return instance;
}

double LJ_Potential::calculateForce(double d) {
    if (d < 0.0) {
        throw std::invalid_argument("Odległość d musi być większa niż 0.");
    }
    else if(d == 0){
        d = 0.01;
    }
    if (d > rcut) {
        return 0.0; // Poza zasięgiem wpływu potencjału Lennarda-Jonesa
    }

    if (forceValues.find(d) != forceValues.end()) {
        return forceValues[d];
    } else {
        double r6 = pow(sigma / d, 6);
        double F = 24.0 * eps * r6 * (2.0 * r6 - 1) / d;

        if (F > 500.0) {
            F = 500.0;
        }
        if (F < -500.0) {
            F = -500.0;
        }

        forceValues[d] = F;

        return F;
    }
}
