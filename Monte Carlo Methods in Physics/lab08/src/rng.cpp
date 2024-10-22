#include "rng.h"

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<double> uniform_dis;

double uniform()
{
    return uniform_dis(gen);
}

double uniform(double a, double b)
{
    std::uniform_real_distribution<double> uniform_dis(a, b);
    return uniform_dis(gen);
}

std::pair<double, double> box_muller(double mu, double s)
{
    double u = uniform();
    double coeff = std::sqrt(-2.0 * std::log(1 - uniform()));

    double x = coeff * std::cos(2.0 * M_PI * u);
    double y = coeff * std::sin(2.0 * M_PI * u);

    return {mu + s * x, mu + s * y};
}