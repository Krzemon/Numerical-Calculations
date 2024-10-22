#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <cstdlib>
#include "generatory.h"

int main() {
    const int N = 1000000; // Liczba punktów
    const int num_bins = 10; // Liczba binów

    double delta_a = 0.5, delta_b = 0.05;
    std::vector<double> zlozony;
    std::vector<double> markova_a;
    std::vector<double> markova_b;
    std::vector<double> elimina;

    markova_a.push_back(uniform());
    markova_b.push_back(uniform());

    for (int i = 0; i < N; ++i) {
        zlozony.push_back(complex_distribution());
        markova_a.push_back(markov_chain(delta_a, markova_a[i]));
        markova_b.push_back(markov_chain(delta_b, markova_b[i]));
        elimina.push_back(elimination_method());
    }

    std::ofstream zlozonyFile("../data/zlozony.dat", std::ios::out);
    std::ofstream markova_aFile("../data/markova_a.dat", std::ios::out);
    std::ofstream markova_bFile("../data/markova_b.dat", std::ios::out);
    std::ofstream eliminaFile("../data/elimina.dat", std::ios::out);
    std::ofstream theoryFile("../data/theory.dat", std::ios::out);

    writeNormalizedHistogramData(zlozonyFile, zlozony, num_bins);  
    writeNormalizedHistogramData(markova_aFile, markova_a, num_bins);
    writeNormalizedHistogramData(markova_bFile, markova_b, num_bins);
    writeNormalizedHistogramData(eliminaFile, elimina, num_bins);

    // Zapisz dane teoretyczne znormalizowane do krzywej gęstości
    for (int i = 0; i <= 100; ++i) {
        double x = static_cast<double>(i) / 100;
        theoryFile << x << " " << f(x) << std::endl;
    }

    zlozonyFile.close();
    markova_aFile.close();
    markova_bFile.close();
    eliminaFile.close();
    theoryFile.close();

    std::vector<double> bin_boundaries(num_bins + 1);
    for (int i = 0; i <= num_bins; ++i) {
        bin_boundaries[i] = static_cast<double>(i) / num_bins;
    }
    
    std::vector<int> n_zlozony(num_bins, 0);
    for (int i = 0; i < N; ++i) {
        int bin_index = static_cast<int>(zlozony[i] * num_bins);
        n_zlozony[bin_index]++;
    }
    double chi2_zlozony = chi2(n_zlozony, bin_boundaries, N);

    std::vector<int> n_markova_a(num_bins, 0);
    for (int i = 0; i < N; ++i) {
        int bin_index = static_cast<int>(markova_a[i] * num_bins);
        n_markova_a[bin_index]++;
    }
    double chi2_markova_a = chi2(n_markova_a, bin_boundaries, N);


    std::vector<int> n_markova_b(num_bins, 0);
    for (int i = 0; i < N; ++i) {
        int bin_index = static_cast<int>(markova_b[i] * num_bins);
        n_markova_b[bin_index]++;
    }
    double chi2_markova_b = chi2(n_markova_b, bin_boundaries, N);

    std::vector<int> n_elimina(num_bins, 0);
    for (int i = 0; i < N; ++i) {
        int bin_index = static_cast<int>(elimina[i] * num_bins);
        n_elimina[bin_index]++;
    }
    double chi2_elimina = chi2(n_elimina, bin_boundaries, N);

    // Próg wartości chi² dla 9 stopni swobody i poziomu istotności 0.05
    double chi2_threshold = 16.9190;

    // Porównanie wartości chi² z wartością graniczną
    std::cout << "Chi2 dla metody złożonej: " << chi2_zlozony << std::endl;
    std::cout << "Chi2 dla metody Łańcucha Markova Delta = 0.5: " << chi2_markova_a << std::endl;
    std::cout << "Chi2 dla metody Łańcucha Markova Delta = 0.05: " << chi2_markova_b << std::endl;
    std::cout << "Chi2 dla metody Eliminacji: " << chi2_elimina << std::endl;
    std::cout << "Wartość graniczna chi2: " << chi2_threshold << std::endl;

    if (chi2_zlozony > chi2_threshold) {
        std::cout << "Hipoteza H0 dla metody złożonej zostaje odrzucona." << std::endl;
    } else {
        std::cout << "Hipoteza H0 dla metody złożonej nie zostaje odrzucona." << std::endl;
    }

    if (chi2_markova_a > chi2_threshold) {
        std::cout << "Hipoteza H0 dla metody Łańcucha Markova Delta = 0.5 zostaje odrzucona." << std::endl;
    } else {
        std::cout << "Hipoteza H0 dla metody Łańcucha Markova Delta = 0.5 nie zostaje odrzucona." << std::endl;
    }

    if (chi2_markova_b > chi2_threshold) {
        std::cout << "Hipoteza H0 dla metody Łańcucha Markova Delta = 0.05 zostaje odrzucona." << std::endl;
    } else {
        std::cout << "Hipoteza H0 dla metody Łańcucha Markova Delta = 0.05 nie zostaje odrzucona." << std::endl;
    }

    if (chi2_elimina > chi2_threshold) {
        std::cout << "Hipoteza H0 dla metody Eliminacji zostaje odrzucona." << std::endl;
    } else {
        std::cout << "Hipoteza H0 dla metody Eliminacji nie zostaje odrzucona." << std::endl;
    }

    return 0;
}
