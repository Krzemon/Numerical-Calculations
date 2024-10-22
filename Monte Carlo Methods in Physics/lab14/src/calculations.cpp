#include "calculations.h"
#include <cmath>
#include <cstdlib>

double uniform()
{
    return (rand() / (double)RAND_MAX);
}

double Psi_T(double r, double a, double c)
{
    return (1 + c * r) * exp(-a * r);
}

double Eps_loc(double r, double a, double c)
{
    return (-a * a * c * r * r + (-a * a + 4 * a * c - 2 * c) * r + 2 * a - 2 * c - 2) / (2 * c * r * r + 2 * r);
}

double Eps_loc_mom(double a, double c, int N, int moment, double delta_r)
{
    double sum = 0.0;
    double r = 0.1;
    for (int i = 0; i < N; i++)
    {
        sum += pow(Eps_loc(r, a, c), moment);
        double U1 = uniform();
        double r_new = r + delta_r * (2 * U1 - 1);
        if (r_new <= 0)
            continue;
        double p_acc = std::min((r_new * r_new * pow(std::abs(Psi_T(r_new, a, c)), 2)) / (r * r * pow(std::abs(Psi_T(r, a, c)), 2)), 1.0);
        double U2 = uniform();
        if (U2 > p_acc)
            continue;
        r = r_new;
    }
    return sum / (double)N;
}
