#ifndef GENERATORY_H
#define GENERATORY_H

#include <vector>
#include <fstream>

double uniform();
double complex_distribution();
double markov_chain(double delta, double x_0);
double elimination_method();
void writeNormalizedHistogramData(std::ofstream& file, const std::vector<double>& data, int num_bins);
double f(double x);
double calculateMaxTheoreticalValue(); // Dodana deklaracja funkcji
double F(double x);
double chi2(const std::vector<int>& n, const std::vector<double>& bins, int N);

#endif // GENERATORY_H
