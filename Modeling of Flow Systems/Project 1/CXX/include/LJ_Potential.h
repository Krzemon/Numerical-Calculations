#ifndef LJ_POTENTIAL_H
#define LJ_POTENTIAL_H

#include "Includes.h"

class LJ_Potential {
public:
    static LJ_Potential& getInstance(double sigma, double eps, double rcut);

    double calculateForce(double d);

private:
    double sigma, eps, rcut;
    std::map<double, double> forceValues;

    LJ_Potential(double sigma, double eps, double rcut);
};

#endif // LJ_POTENTIAL_H
