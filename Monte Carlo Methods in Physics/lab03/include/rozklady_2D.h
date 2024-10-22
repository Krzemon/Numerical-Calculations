#ifndef GENERATORY_H
#define GENERATORY_H

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <fstream>

const int n = 1e4;

inline double uniform() {
    return static_cast<double>(std::rand()) / RAND_MAX;
}

void Box_Muller(std::vector<double> &X, std::vector<double> &Y);

void Circle_Norm(std::vector<double> &X, std::vector<double> &Y);

void Variance(std::vector<double> &X, std::vector<double> &Y, double cov[2][2]);

void Print_Cov(double cov[2][2]);

void Transformance(double scale1, double scale2, double angle, double vec1[2], double vec2[2],
                   std::vector<double> &X, std::vector<double> &Y,
                   std::vector<double> &X_result, std::vector<double> &Y_result);
                

#endif // GENERATORY_H
