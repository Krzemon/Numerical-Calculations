#ifndef RNG_H
#define RNG_H

#include <random>
#include <tuple>

extern std::random_device rd;
extern std::mt19937 gen;
extern std::uniform_real_distribution<double> uniform_dis;

double uniform();
double uniform(double a, double b);

std::pair<double,double> box_muller(double mu=0.0, double s=1.0);

#endif