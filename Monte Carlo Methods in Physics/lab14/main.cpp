#include "calculations.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

// Define constants
const int N = 1000000; // 10^6
const double delta_r = 0.1;
const double a_min = 0.3, a_max = 1.2;
const double c_min = -0.7, c_max = 0.3;
const double delta_a = 0.02, delta_c = 0.02;
const int M = 200;
const double r_max = 8.0;

int main()
{
    int n_a = static_cast<int>((a_max - a_min) / delta_a) + 1;
    int n_c = static_cast<int>((c_max - c_min) / delta_c) + 1;

    std::vector<std::vector<double>> eps_mom1(n_a, std::vector<double>(n_c));
    std::vector<std::vector<double>> eps_mom2(n_a, std::vector<double>(n_c));
    std::vector<std::vector<double>> sigma(n_a, std::vector<double>(n_c));
    std::vector<std::vector<double>> log_sigma(n_a, std::vector<double>(n_c));

    for (int a_idx = 0; a_idx < n_a; a_idx++)
    {
        for (int c_idx = 0; c_idx < n_c; c_idx++)
        {
            double a = a_min + a_idx * delta_a;
            double c = c_min + c_idx * delta_c;
            eps_mom1[a_idx][c_idx] = Eps_loc_mom(a, c, N, 1, delta_r);
            eps_mom2[a_idx][c_idx] = Eps_loc_mom(a, c, N, 2, delta_r);
            sigma[a_idx][c_idx] = sqrt(eps_mom2[a_idx][c_idx] - pow(eps_mom1[a_idx][c_idx], 2));
            log_sigma[a_idx][c_idx] = log(sigma[a_idx][c_idx]);
        }
    }

    // Output results to files
    std::ofstream eps("../data/eps_mom1.dat");
    std::ofstream sigma_file("../data/sigma.dat");
    std::ofstream log_sigma_file("../data/log_sigma.dat");

    if (eps.is_open() && sigma_file.is_open() && log_sigma_file.is_open()){
        for (int a_idx = 0; a_idx < n_a; a_idx++)
        {
            for (int c_idx = 0; c_idx < n_c; c_idx++)
            {
                eps << eps_mom1[a_idx][c_idx] << " ";
                sigma_file << sigma[a_idx][c_idx] << " ";
                log_sigma_file << log_sigma[a_idx][c_idx] << " ";
            }
            eps << std::endl;
            sigma_file << std::endl;
            log_sigma_file << std::endl;
        }
        eps.close();
        sigma_file.close();
        log_sigma_file.close();
        std::cout << "Data saved" << std::endl;
    }
    else
    {
        std::cerr << "Unable to open file for writing" << std::endl;
        return 1;
    }

    // Generate histogram for ground state (a = 1, c = 0)
    double a = 1.0;
    double c = 0.0;
    std::vector<double> dist(M, 0.0);

    double r = 0.1;
    for (int i = 0; i < N; i++)
    {
        double U1 = uniform();
        double r_new = r + delta_r * (2 * U1 - 1);
        if (r_new <= 0)
            continue;
        double p_acc = std::min((r_new * r_new * pow(std::abs(Psi_T(r_new, a, c)), 2)) / (r * r * pow(std::abs(Psi_T(r, a, c)), 2)), 1.0);
        double U2 = uniform();
        if (U2 > p_acc)
            continue;
        r = r_new;

        if (r <= r_max)
        {
            int k = static_cast<int>(floor(r / (r_max / M)));
            dist[k] += 1.0 / (N * (r_max / M));
        }
    }

    // Save histogram data
    std::ofstream hist_file("../data/histogram.dat");
    if (hist_file.is_open())
    {
        for (int k = 0; k < M; k++)
        {
            double r_bin = k * (r_max / M);
            hist_file << r_bin << " " << dist[k] << std::endl;
        }
        hist_file.close();
        std::cout << "Histogram data saved" << std::endl;
    }
    else
    {
        std::cerr << "Unable to open histogram file for writing" << std::endl;
        return 1;
    }

    return 0;
}