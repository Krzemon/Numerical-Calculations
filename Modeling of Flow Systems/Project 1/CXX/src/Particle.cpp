#include "Particle.h"

Particle::Particle(double x, double y, double m) : xpos(x), ypos(y), mass(m), xvel(0.0), yvel(0.0), xfor(0.0), yfor(0.0) {}

Particle::Particle(double x, double y, double m, int xD) : xpos(x), ypos(y), mass(m), xvel(0.0), yvel(0.0), xfor(0.0), yfor(0.0) {
    // Assuming xD is a flag to set the initial velocity in the x-direction
    if (xD == 1) {
        xvel = 1.0;  // You can set the desired initial velocity
    }
}

Particle::Particle(double x, double y, double m, double vx, double vy) : xpos(x), ypos(y), mass(m), xvel(vx), yvel(vy), xfor(0.0), yfor(0.0) {}

void Particle::integrate(double dt, double Fc, double ixh) {
    // Assuming dt is the time step, Fc is a force coefficient, and ixh is some parameter
    // Update velocity using Verlet integration
    xvel_1 = xvel + 0.5 * xfor / mass * dt;
    yvel_1 = yvel + 0.5 * yfor / mass * dt;

    // Update position
    xpos += xvel_1 * dt + 0.5 * xfor / mass * dt * dt;
    ypos += yvel_1 * dt + 0.5 * yfor / mass * dt * dt;

    // Update forces (assuming Fc and ixh are used in the force calculation)
    xfor = -Fc * xvel - ixh * xvel_1;
    yfor = -Fc * yvel - ixh * yvel_1;

    // Update velocity again
    xvel_2 = xvel_1 + 0.5 * xfor / mass * dt;
    yvel_2 = yvel_1 + 0.5 * yfor / mass * dt;

    xvel = xvel_2;
    yvel = yvel_2;
}