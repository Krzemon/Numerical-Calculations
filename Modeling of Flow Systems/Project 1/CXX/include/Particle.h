#ifndef PARTICLE_H
#define PARTICLE_H

#include "Includes.h"

class Particle {
public:
    double xpos, ypos, xvel, yvel, xvel_1, yvel_1, xvel_2, yvel_2, xfor, yfor, mass;

    Particle(double x, double y, double m);
    Particle(double x, double y, double m, int xD);
    Particle(double x, double y, double m, double vx, double vy);
    void integrate(double dt, double Fc, double ixh);
};

#endif // PARTICLE_H