#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <cstdlib>
#include "generatory.h"

static double g1 = 4.0 / 5, g2 = 1.0 / 5;

double uniform() {
    return static_cast<double>(std::rand()) / RAND_MAX;
}

double complex_distribution() {
    double u1 = uniform();
    double u2 = uniform();

    if (u1 <= g1)
        return u1;
    else
        return std::sqrt(1 - std::sqrt(1 - u2));
}

double markov_chain(double delta, double x_0) {
    double x_new = x_0 + delta * (2 * uniform() - 1);
    double p_acc = std::min(1.0, f(x_new) / f(x_0));

    if (uniform() <= p_acc && x_new >= 0 && x_new <= 1)
        return x_new;
    else
        return x_0;
}

double elimination_method() {
    double x = uniform();
    double y = 1.15 * uniform();
    double p = f(x);

    if (y <= p)
        return x;
    else
        return elimination_method();
}

void writeNormalizedHistogramData(std::ofstream& file, const std::vector<double>& data, int num_bins) {
    std::vector<int> histogram(num_bins, 0);

    // Oblicz histogram
    for (auto& value : data) {
        int bin_index = static_cast<int>(value * num_bins);
        bin_index = std::min(std::max(bin_index, 0), num_bins - 1); // Ogranicz indeks binu
        histogram[bin_index]++;
    }

    int max_count = *std::max_element(histogram.begin(), histogram.end());

    for (int i = 0; i < num_bins; ++i) {
        double bin_start = static_cast<double>(i) / num_bins;
        double bin_end = static_cast<double>(i + 1) / num_bins;
        double bin_center = (bin_start + bin_end) / 2.0 + 1.0 / (2 * num_bins);

        // Znormalizowana wysokość słupka
        double normalized_height = static_cast<double>(histogram[i]) / max_count;
        double max_theoretical = calculateMaxTheoreticalValue();
        // normowanie do wysokoścu wykresu teoretycznego
        file << bin_center << " " << normalized_height * max_theoretical<< std::endl;
    }
}

double f(double x) {
    return 4.0 / 5 * (1 + x - std::pow(x, 3.0));
}

double F(double x) {
    return 4.0 / 5 * (x + std::pow(x, 2.0) / 2 - std::pow(x, 4.0) / 4);
}

double calculateMaxTheoreticalValue() {
    std::vector<double> theory_values;
    for (int i = 0; i <= 100; ++i) {  // Zapisz więcej punktów dla krzywej teoretycznej
        double x = static_cast<double>(i) / 100;
        theory_values.push_back(f(x));
    }
    return *std::max_element(theory_values.begin(), theory_values.end());
}

double chi2(const std::vector<int>& n, const std::vector<double>& bins, int N) {
    double chi2_value = 0.0;
    for (size_t i = 0; i < n.size(); ++i) {
        double expected = (F(bins[i + 1]) - F(bins[i])) * N;
        chi2_value += std::pow(n[i] - expected, 2) / expected;
    }
    return chi2_value;
}
