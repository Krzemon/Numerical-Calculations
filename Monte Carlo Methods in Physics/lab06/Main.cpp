#include "wiener.h"
#include "particle_translation.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <random>    // dla generatora rozkładu normalnego
#include <string>
#include <cstdlib>

struct Particle {
    double x, y;
    int theta;
    double l;
};

void gaus01(double& x, double& y) {
    double u1, u2;
    u1 = uniform();
    u2 = uniform();
    x = sqrt(-2.0 * log(1.0 - u1)) * cos(2.0 * M_PI * u2);
    y = sqrt(-2.0 * log(1.0 - u1)) * sin(2.0 * M_PI * u2);
}


struct ParticleSystem {
    int N;
    double R, xr, yr, ra, xa, ya, D, sigmaDeltaT;
    bool closed;
    Particle* particles;
    
    ParticleSystem(int Ne, double xs, double ys, double Re, double xre, double yre, double rae, double xae, double yae, double De) {
        N = Ne;
        R = Re;
        xr = xre;
        yr = yre;
        ra = rae;
        xa = xae;
        ya = yae;
        D = De;
        closed = true;
        particles = new Particle[N];
        for (int i = 0; i < N; i++) {
            particles[i].x = xs;
            particles[i].y = ys;
            particles[i].theta = 0;
        }
    }

    void loop(double tmax, double deltaT, double omega, std::ofstream& activeFile, double xa, double ya, double Ra) {

        int M = int(round(tmax / deltaT));
        double deltaN = omega * deltaT;
        sigmaDeltaT = sqrt(2.0 * D * deltaT);
        int counter;
        double* xPositions = new double[N];
        double* yPositions = new double[N];

        for (int j = 0; j < M; j++) {
            int activeCount = 0;
            for (int i = 0; i < N; i++) {
                if (particles[i].theta == 1) {
                    activeCount++;
                }
            }



            activeFile << j * deltaT << " " << activeCount << std::endl;


            for (int i = 0; i < N; i++) {
                xPositions[i] = particles[i].x;
                yPositions[i] = particles[i].y;
            }
            


            for (int i = 0; i < N; i++) {
                counter = deltaN;
                if (particles[i].theta == 1) {
                    particles[i].l = 1.0;
                    double deltaX, deltaY;
                    gaus01(deltaX, deltaY);
                    double x2 = particles[i].x + deltaX * sigmaDeltaT, y2 = particles[i].y + deltaY * sigmaDeltaT;
                    do {
                        particle_translation(particles[i].x, particles[i].y, x2, y2, xr, yr, R, xa, ya, ra, particles[i].theta, particles[i].l);
                    } while (particles[i].l > 1E-6);
                }
                // dopływ cząstek
                if(counter > 0 && particles[i].theta == 0){
                    particles[i].theta = 1;
                    counter--;
                }   
            }



            // absorbent pochłania 
            for (int i = 0; i < N; i++) {
                if (xPositions[i] > xa - Ra && xPositions[i] < xa + Ra && yPositions[i] > ya - Ra && yPositions[i] < ya + Ra && particles[i].theta == 1) {
                    particles[i].theta = 0;
                    activeCount--;
                }
            }

            if (j == 1 || j == 10 || j == 100 || j == 1000) {
                std::ofstream file("../data/position_" + std::to_string(j) + ".dat", std::ios::out | std::ios::trunc);
                for (int i = 0; i < N; ++i) {
                    if (particles[i].theta == 1 && (xPositions[i] < xa - Ra || xPositions[i] > xa + Ra || yPositions[i] < ya - Ra || yPositions[i] > ya + Ra )) {
                        file << xPositions[i] << " " << yPositions[i] << std::endl;
                    }
                }
                file.close();
            }

        }

        delete[] xPositions;
        delete[] yPositions;
    }
    };



int main() {

    // zadanie 1
    // wiener();

    // zadaie 2

    double tmax = t_max;
    double omega_values[3] = {10.0, 50.0, 100.0};
    double xr = 0.0, yr = 0.0, R = 5.0;
    double xz = -4.5, yz = 0.0;
    double xa = 3.0, ya = 0.0;
    double Ra_values[2] = {0.1, 0.5};

    double omegaValues[3] = {10, 50, 100};
    double RaValues[2] = {0.1, 0.5};

    std::ofstream active("../data/active.dat", std::ios::out | std::ios::trunc);
    ParticleSystem system = ParticleSystem(10000, xz, yz, R, xr, yr, RaValues[1], xa, ya, 1.0);
    system.loop(1000, 0.1, omegaValues[1], active, xa, ya, Ra_values[1]);
    active.close();


    return 0;
}