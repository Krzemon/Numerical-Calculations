#ifndef Fullerene_h
#define Fullerene_h

#include <vector>
#include <tuple>
#include <string_view>
#include <cmath>

#include <Eigen/Dense>

using pos = Eigen::Vector3d;
using pos_list = Eigen::VectorX<pos>;

double calc_r(const pos &first, const pos &second);
double calc_r(size_t i, size_t j, const pos_list &atoms);

pos to_cartesian(const pos &spherical);
pos to_spherical(const pos &cartesian);

class Fullerene
{
public:
    Fullerene(size_t size, double r);
    Fullerene(pos_list &&cartesian, pos_list &&spherical);
    static Fullerene from_file(std::string_view filename);

    void restrict_quadruple_bonds(double xi = 10.0);
    void allow_quadruple_bonds();

    void write_positions(std::string_view filename) const;

    size_t get_size() const
    {
        return atoms.size();
    }

    double get_energy() const
    {
        return V_tot;
    }

    double get_r_avg() const
    {
        return r_avg;
    }

    Eigen::VectorXd calc_PCF(size_t M, double rmax_coeff = 2.5);

    void try_shifting(double beta, double w_r, double w_phi, double w_theta, double W_all, bool use_full_potential = false);

private:
    struct BrennerPotential
    {
        BrennerPotential() = delete;
        BrennerPotential(size_t size,
                         bool limit_bonds = false,
                         double xi = 10.0,
                         double R0 = 1.315,
                         double R1 = 1.70,
                         double R2 = 2.00,
                         double De = 6.325,
                         double S = 1.29,
                         double lambda = 1.5,
                         double delta = 0.80469,
                         double a0 = 0.011304,
                         double c0 = 19,
                         double d0 = 2.5);

        double operator()(const pos_list &atoms);

        size_t n;
        bool limit_bonds;
        double xi;

        double R0; // angstrom
        double R1; // angstrom
        double R2; // angstrom
        double De; // eV
        double S;
        double lambda; // 1/angstrom
        double delta;
        double a0;
        double c0;
        double d0;

        double V_coeff;
        double g_prep;


        Eigen::MatrixXd r;
        Eigen::MatrixXd f;
        Eigen::MatrixXd V_R;
        Eigen::MatrixXd V_A;
        Eigen::MatrixXd B_avg;
        Eigen::VectorXd V;

        void fill_r(const pos_list &atoms);
        void fill_f_VR_VA();
        void fill_B_avg(const pos_list &atoms);
        void fill_V();

        double calc_f(double r_ij) const
        {
            if (r_ij <= R1)
                return 1.0;

            if (r_ij <= R2)
                return 0.5 * (1.0 + std::cos(M_PI * (r_ij - R1) / (R2 - R1)));

            return 0.0;
        }

        double calc_V_R(double r_ij) const
        {
            return De / (S - 1.0) * std::exp(-std::sqrt(2.0 * S) * lambda * (r_ij - R0));
        }

        double calc_V_A(double r_ij) const
        {
            return De / (S - 1.0) * S * std::exp(-std::sqrt(2.0 / S) * lambda * (r_ij - R0));
        }

        double calc_B_avg(size_t i, size_t j, const pos_list &atoms) const;
        double calc_V(size_t i, const pos_list &atoms) const;

    } BP;

    pos_list atoms; //cartesian coordinates
    pos_list atoms_spherical; //spherical coordinates
    bool limit_bonds;

    double V_tot{};
    double r_avg{};

    void calc_r_avg();
};

#endif