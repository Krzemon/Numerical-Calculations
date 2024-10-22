#include "Fullerene.h"

#include <iostream>
#include <fstream>
#include <cmath>

int main()
{
    bool zad2 = 1;
    bool zad3 = 1;
    bool zad4 = 1;
    bool zad5 = 1;
    bool zad6 = 1;
    bool zad7 = 1;

    size_t n = 60;
    double beta_min = 1.0;
    double beta_max = 100.0;
    double p = 2.0;
    double w_r = 1e-4;
    double w_phi = 0.05;
    double w_theta = 0.05;
    double W_all = 1e-4;

    double nn_coeff = 1.0;

    size_t N = 2e5;

    double r_init = 3.5;
    size_t M = 100;

    // zad2
    if (zad2)
    {
        Fullerene C60 = Fullerene::from_file("../C60.dat");

        std::cout << "zad1: " << std::endl;
        std::cout << "V_tot = " << C60.get_energy() << " eV" << std::endl;
        std::cout << "E_b = " << C60.get_energy() / C60.get_size() << " eV" << std::endl;
        std::cout << "r_avg = " << C60.get_r_avg() << " A" << std::endl;

        C60.write_positions("data/zad2_C60.dat");
        std::ofstream PCF_file("data/zad2_PCF.dat");

        PCF_file << C60.calc_PCF(M);
    }

    // zad3
    if (zad3)
    {
        Fullerene C60(n, r_init);

        std::ofstream zad3_file("data/zad3.dat");

        for (size_t it = 1; it <= N; ++it)
        {
            double beta = beta_min + (beta_max - beta_min) * std::pow(static_cast<double>(it) / N, p);

            C60.try_shifting(beta, w_r, w_phi, w_theta, W_all);

            if (it % 100 == 0)
            {
                zad3_file << it << " " << beta << " " << C60.get_energy() << " " << C60.get_r_avg() << std::endl;
            }
        }

        C60.write_positions("data/zad3_C60.dat");
        std::ofstream PCF_file("data/zad3_PCF.dat");

        PCF_file << C60.calc_PCF(M);
    }

    // zad4
    if (zad4)
    {

        Fullerene C60(n, r_init);
        C60.restrict_quadruple_bonds();

        std::ofstream zad4_file("data/zad4.dat");

        for (size_t it = 1; it <= N; ++it)
        {
            double beta = beta_min + (beta_max - beta_min) * std::pow(static_cast<double>(it) / N, p);

            C60.try_shifting(beta, w_r, w_phi, w_theta, W_all);

            if (it % 100 == 0)
            {
                zad4_file << it << " " << beta << " " << C60.get_energy() << " " << C60.get_r_avg() << std::endl;
            }
        }

        C60.write_positions("data/zad4_C60.dat");
        std::ofstream PCF_file("data/zad4_PCF.dat");

        PCF_file << C60.calc_PCF(M);
    }

    // zad5
    if (zad5)
    {
        Fullerene C60(n, 2.5);
        C60.restrict_quadruple_bonds();

        std::ofstream zad5_file("data/zad5.dat");

        for (size_t it = 1; it <= N; ++it)
        {
            double beta = beta_min + (beta_max - beta_min) * std::pow(static_cast<double>(it) / N, p);

            C60.try_shifting(beta, w_r, w_phi, w_theta, W_all);

            if (it % 100 == 0)
            {
                zad5_file << it << " " << beta << " " << C60.get_energy() << " " << C60.get_r_avg() << std::endl;
            }
        }

        C60.write_positions("data/zad5_C60.dat");
        std::ofstream PCF_file("data/zad5_PCF.dat");

        PCF_file << C60.calc_PCF(M);
    }

    if (zad6)
    {
        // beta_min, beta_max, p, w_r, w_phi, w_theta, W_all
        std::vector<std::tuple<double, double, double, double, double, double>> params = {
            {beta_min / 10, beta_max, p, w_r, w_phi, w_theta},
            {beta_min, beta_max * 5, p, w_r, w_phi, w_theta},
            {beta_min, beta_max, 3.0, w_r, w_phi, w_theta},
            {beta_min, beta_max, p, w_r * 5, w_phi, w_theta},
            {beta_min, beta_max, p, w_r, w_phi * 5, w_theta},
            {beta_min, beta_max, p, w_r, w_phi, w_theta * 5}};

        for (int i = 0; i < params.size(); i++)
        {
            Fullerene C60(n, 2.5);
            C60.restrict_quadruple_bonds();

            std::ofstream zad6_file("data/zad6_" + std::to_string(i) + ".dat");

            auto [beta_min, beta_max, p, w_r, w_phi, w_theta] = params[i];

            for (size_t it = 1; it <= N; ++it)
            {
                double beta = beta_min + (beta_max - beta_min) * std::pow(static_cast<double>(it) / N, p);

                C60.try_shifting(beta, w_r, w_phi, w_theta, W_all);

                if (it % 100 == 0)
                {
                    zad6_file << it << " " << beta << " " << C60.get_energy() << " " << C60.get_r_avg() << std::endl;
                }
            }

            C60.write_positions("data/zad6_C60_" + std::to_string(i) + ".dat");
            std::ofstream PCF_file("data/zad6_PCF_" + std::to_string(i) + ".dat");

            PCF_file << C60.calc_PCF(M);
        }
    }

    if (zad7)
    {
        std::ofstream zad7_file("data/zad7.dat");

        int n_min = 20;
        int n_max = 60;

        double r_init_min = 1.3;
        double r_init_max = 2.5;

        for (int n = n_min; n <= n_max; n++)
        {
            Fullerene Cn(n, r_init_min + (r_init_max - r_init_min) * static_cast<double>(n - n_min) / (n_max - n_min));
            Cn.restrict_quadruple_bonds();

            double E_sum = 0.0;
            double E2_sum = 0.0;
            double r_sum = 0.0;
            double r2_sum = 0.0;
            size_t n_sum = 100;

            double beta_max = 200.0;

            for (size_t it = 1; it <= N; ++it)
            {
                double beta = beta_min + (beta_max - beta_min) * std::pow(static_cast<double>(it) / N, p);
                Cn.try_shifting(beta, w_r, w_phi, w_theta, W_all);

                if (it > N - n_sum)
                {
                    double E = Cn.get_energy() / n;
                    E_sum += E;
                    E2_sum += E * E;

                    double r = Cn.get_r_avg();
                    r_sum += r;
                    r2_sum += r * r;
                }
            }

            double E_avg = E_sum / n_sum;
            double E2_avg = E2_sum / n_sum;
            double r_avg = r_sum / n_sum;
            double r2_avg = r2_sum / n_sum;

            double E_std = std::sqrt(E2_avg - E_avg * E_avg);
            double r_std = std::sqrt(r2_avg - r_avg * r_avg);

            zad7_file << n << " " << E_avg << " " << E_std << " " << r_avg << " " << r_std << std::endl;
        }
    }
}
