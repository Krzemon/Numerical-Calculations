#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <iomanip>

// Funkcja generująca losową liczbę rzeczywistą z zakresu [0, 1)
double uniform();

// Funkcja generująca punkty na kole z rozkładem jednorodnym
void sferycznie_konturowany(std::vector<double> &X, std::vector<double> &Y, double R);

// Funkcja przekształcająca punkty na obwód okregu + przesunięcie
void circle_Norm(std::vector<double> &X, std::vector<double> &Y);

// Funkcja obliczająca wspólną powierzchnię przekrywających się kół o promieniach R_alpha i R_beta
double calculateOverlapArea(const std::vector<double> &X, const std::vector<double> &Y, double R_alpha, double R_beta, double x_A, int n);

// Funkcja przesuwająca środek okręgu o wektor (x_A, y_A)
void translation(std::vector<double> &X, std::vector<double> &Y, double x_A, double y_A);


void drawPoints(const std::vector<double> &X, const std::vector<double> &Y, const std::string &filename);
double calculateVariance(double mean, double secondMoment, int N);
double calculateStandardDeviation(double mean, double secondMoment, int N);
double calculateCircleArea(double R_alpha);


#endif // GEOMETRY_H