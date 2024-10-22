#include "rng.h"
#include "Fullerene.h"

#include <fstream>
#include <iostream>

double calc_r(const pos &first, const pos &second)
{
    return (first - second).norm();
}

double calc_r(size_t i, size_t j, const pos_list &atoms)
{
    return calc_r(atoms(j), atoms(i));
}

pos to_cartesian(const pos &spherical)
{
    double r = spherical.x();
    double phi = spherical.y();
    double theta = spherical.z();

    double x = r * std::sin(theta) * std::cos(phi);
    double y = r * std::sin(theta) * std::sin(phi);
    double z = r * std::cos(theta);

    return pos(x, y, z);
}

pos to_spherical(const pos &cartesian)
{
    double x = cartesian.x();
    double y = cartesian.y();
    double z = cartesian.z();

    double r = std::sqrt(x * x + y * y + z * z);
    double phi = std::atan2(y, x);
    double theta = std::acos(z / r);

    return pos(r, phi, theta);
}

Fullerene::BrennerPotential::BrennerPotential(size_t size,
                                              bool limit_bonds,
                                              double xi,
                                              double R0,
                                              double R1,
                                              double R2,
                                              double De,
                                              double S,
                                              double lambda,
                                              double delta,
                                              double a0,
                                              double c0,
                                              double d0) : n(size),
                                                           limit_bonds(limit_bonds),
                                                           xi(xi),
                                                           R0(R0),
                                                           R1(R1),
                                                           R2(R2),
                                                           De(De),
                                                           S(S),
                                                           lambda(lambda),
                                                           delta(delta),
                                                           a0(a0),
                                                           c0(c0),
                                                           d0(d0),
                                                           V_coeff(De / (S - 1.0)),
                                                           g_prep(1 + c0 * c0 / d0 / d0),
                                                           r(Eigen::MatrixXd::Zero(size, size)),
                                                           f(Eigen::MatrixXd::Zero(size, size)),
                                                           V_R(Eigen::MatrixXd::Zero(size, size)),
                                                           V_A(Eigen::MatrixXd::Zero(size, size)),
                                                           B_avg(Eigen::MatrixXd::Zero(size, size)),
                                                           V(Eigen::VectorXd::Zero(size))
{
}

double Fullerene::BrennerPotential::operator()(const pos_list &atoms)
{
    fill_r(atoms);
    fill_f_VR_VA();
    fill_B_avg(atoms);
    fill_V();

    return 0.5 * V.sum();
}

void Fullerene::BrennerPotential::fill_r(const pos_list &atoms)
{
    for (size_t i = 0; i < n; ++i)
        for (size_t j = i + 1; j < n; ++j)
            r(i, j) = r(j, i) = calc_r(i, j, atoms);
}

void Fullerene::BrennerPotential::fill_f_VR_VA()
{
    double r_ij = 0.0;
    double arg = 0.0;

    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            r_ij = r(i, j);

            if (r_ij > R2 || i == j)
            {
                f(i, j) = f(j, i) = 0.0;
                V_R(i, j) = V_R(j, i) = 0.0;
                V_A(i, j) = V_A(j, i) = 0.0;

                continue;
            }

            if (r_ij <= R1)
            {
                f(i, j) = f(j, i) = 1.0;
            }
            else
            {
                f(i, j) = f(j, i) = 0.5 * (1.0 + std::cos(M_PI * (r_ij - R1) / (R2 - R1)));
            }

            arg = lambda * (r_ij - R0);

            V_R(i, j) = V_R(j, i) = V_coeff * std::exp(-std::sqrt(2.0 * S) * arg);
            V_A(i, j) = V_A(j, i) = V_coeff * S * std::exp(-std::sqrt(2.0 / S) * arg);
        }
    }
}

void Fullerene::BrennerPotential::fill_B_avg(const pos_list &atoms)
{
    double r_ij = 0.0;

    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = i + 1; j < n; ++j)
        {
            r_ij = r(i, j);

            if (r_ij > R2 || i == j)
            {
                B_avg(i, j) = B_avg(j, i) = 0.0;
                continue;
            }

            pos vr_ij = atoms(j) - atoms(i);
            pos vr_ji = -vr_ij;

            double xi_ij = 0.0;
            double xi_ji = 0.0;

            for (size_t k = 0; k < n; ++k)
            {
                double r_ik = r(i, k);

                if (r_ik > R2 || k == i || k == j)
                    continue;

                double cos_theta_ijk = vr_ij.dot(atoms(k) - atoms(i)) / r_ij / r_ik;

                if (limit_bonds && cos_theta_ijk > 0.0)
                {
                    xi_ij = xi;
                    break;
                }

                xi_ij += f(i, k) * a0 * (g_prep - c0 * c0 / (d0 * d0 + (1.0 + cos_theta_ijk) * (1.0 + cos_theta_ijk)));
            }

            for (size_t k = 0; k < n; ++k)
            {
                double r_jk = r(j, k);

                if (r_jk > R2 || k == i || k == j)
                    continue;

                double cos_theta_jik = vr_ji.dot(atoms(k) - atoms(j)) / r_ij / r_jk;

                if (limit_bonds && cos_theta_jik > 0.0)
                {
                    // xi_ji = xi;
                    // break;

                    xi_ji += xi;
                }
                else
                {
                    xi_ji += f(j, k) * a0 * (g_prep - c0 * c0 / (d0 * d0 + (1.0 + cos_theta_jik) * (1.0 + cos_theta_jik)));
                }
            }

            B_avg(i, j) = B_avg(j, i) = 0.5 * (std::pow(1.0 + xi_ij, -delta) + std::pow(1.0 + xi_ji, -delta));
        }
    }
}

void Fullerene::BrennerPotential::fill_V()
{
    V.setZero();

    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            if (r(i, j) > R2 || i == j)
                continue;

            V(i) += f(i, j) * (V_R(i, j) - B_avg(i, j) * V_A(i, j));
        }
    }
}

double Fullerene::BrennerPotential::calc_V(size_t i, const pos_list &atoms) const
{
    double V = 0.0;
    double r_ij = 0.0;

    for (size_t j = 0; j < n; ++j)
    {
        r_ij = calc_r(atoms(i), atoms(j));

        if (r_ij > R2 || i == j)
            continue;

        V += calc_f(r_ij) * (calc_V_R(r_ij) - calc_B_avg(i, j, atoms) * calc_V_A(r_ij));
    }

    return V;
}

double Fullerene::BrennerPotential::calc_B_avg(size_t i, size_t j, const pos_list &atoms) const
{
    double r_ij = calc_r(atoms(i), atoms(j));

    if (r_ij > R2 || i == j)
    {
        return 0.0;
    }

    pos vr_ij = atoms(j) - atoms(i);
    pos vr_ji = -vr_ij;

    double xi_ij = 0.0;
    double xi_ji = 0.0;

    for (size_t k = 0; k < n; ++k)
    {
        double r_ik = calc_r(atoms(i), atoms(k));

        if (r_ik > R2 || k == i || k == j)
            continue;

        double cos_theta_ijk = vr_ij.dot(atoms(k) - atoms(i)) / r_ij / r_ik;

        if (limit_bonds && cos_theta_ijk > 0.0)
        {
            xi_ij = xi;
            break;
        }

        xi_ij += calc_f(r_ij) * a0 * (g_prep - c0 * c0 / (d0 * d0 + (1.0 + cos_theta_ijk) * (1.0 + cos_theta_ijk)));
    }

    for (size_t k = 0; k < n; ++k)
    {
        double r_jk = calc_r(atoms(j), atoms(k));

        if (r_jk > R2 || k == i || k == j)
            continue;

        double cos_theta_jik = vr_ji.dot(atoms(k) - atoms(j)) / r_ij / r_jk;

        if (limit_bonds && cos_theta_jik > 0.0)
        {
            // xi_ji = xi;
            // break;

            xi_ji += xi;
        }
        else
        {
            xi_ji += calc_f(r_jk) * a0 * (g_prep - c0 * c0 / (d0 * d0 + (1.0 + cos_theta_jik) * (1.0 + cos_theta_jik)));
        }
    }

    return 0.5 * (std::pow(1.0 + xi_ij, -delta) + std::pow(1.0 + xi_ji, -delta));
}

Fullerene::Fullerene(size_t size, double r) : BP(size), atoms(size), atoms_spherical(size), r_avg(r)
{
    for (size_t i = 0; i < size; i++)
    {
        double phi = uniform() * 2.0 * M_PI;
        double theta = uniform() * M_PI;

        pos atom_spherical(r, phi, theta);

        atoms(i) = to_cartesian(atom_spherical);
        atoms_spherical(i) = std::move(atom_spherical);
    }

    V_tot = BP(atoms);
}

Fullerene::Fullerene(pos_list &&cartesian, pos_list &&spherical) : BP(cartesian.size()), atoms(cartesian), atoms_spherical(spherical)
{
    V_tot = BP(atoms);
    calc_r_avg();
}

Fullerene Fullerene::from_file(std::string_view filename)
{
    std::ifstream file(filename.data());
    std::string token;

    std::vector<double> coordinates;
    coordinates.reserve(3 * 60);

    while (file >> token)
    {
        coordinates.push_back(std::stod(token));
    }

    pos_list cartesian(coordinates.size() / 3);
    pos_list spherical(coordinates.size() / 3);

    for (long int i = 0; i < cartesian.size(); ++i)
    {
        pos atom(coordinates[3 * i], coordinates[3 * i + 1], coordinates[3 * i + 2]);

        spherical(i) = to_spherical(atom);
        cartesian(i) = std::move(atom);
    }

    return Fullerene(std::move(cartesian), std::move(spherical));
}

void Fullerene::restrict_quadruple_bonds(double xi)
{
    BP.limit_bonds = true;
    BP.xi = xi;
}

void Fullerene::allow_quadruple_bonds()
{
    BP.limit_bonds = false;
}

void Fullerene::write_positions(std::string_view filename) const
{
    std::ofstream file(filename.data());

    for(pos atom : atoms)
    {
        file << atom.x() << " " << atom.y() << " " << atom.z() << std::endl;
    }
}

Eigen::VectorXd Fullerene::calc_PCF(size_t M, double rmax_coeff)
{
    double r_max = rmax_coeff * r_avg;
    double dr = r_max / M;

    Eigen::VectorXd PCF = Eigen::VectorXd::Zero(M);

    size_t n = atoms.size();

    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = i + 1; j < n; ++j)
        {
            double r_ij = BP.r(i, j);
            size_t m = r_ij / dr;

            if (m < M)
                PCF[m] += 2.0 * 4.0 * M_PI * r_avg * r_avg / (n * n * 2.0 * M_PI * r_ij * dr);
        }
    }

    return PCF;
}

void Fullerene::try_shifting(double beta, double w_r, double w_phi, double w_theta, double W_all, bool use_full_potential)
{
    for (long int i = 0; i < atoms.size(); ++i)
    {
        pos atom_old = atoms(i);
        pos atom_old_spherical = atoms_spherical(i);

        double r_new = atom_old_spherical.x() * (1.0 + (2.0 * uniform() - 1) * w_r);
        double phi_new = atom_old_spherical.y() * (1.0 + (2.0 * uniform() - 1) * w_phi);
        double theta_new = atom_old_spherical.z() * (1.0 + (2.0 * uniform() - 1) * w_theta);

        if (phi_new < 0.0)
            phi_new += 2 * M_PI;

        if (phi_new > 2.0 * M_PI)
            phi_new -= 2.0 * M_PI;

        if (theta_new < 0.0 || theta_new > M_PI)
            theta_new = atom_old_spherical.z();

        pos atom_new_spherical(r_new, phi_new, theta_new);

        atoms(i) = std::move(to_cartesian(atom_new_spherical));

        double dV = 0.0;

        if (use_full_potential)
        {
            dV = BP(atoms) - V_tot;
        }
        else
        {
            dV = BP.calc_V(i, atoms) - BP.V(i);
        }

        if (dV <= 0.0)
        {
            atoms_spherical(i) = std::move(atom_new_spherical);
            continue;
        }

        double p = std::exp(-beta * dV);

        if (uniform() <= p)
        {
            atoms_spherical(i) = std::move(atom_new_spherical);
        }
        else
        {
            atoms(i) = std::move(atom_old);
        }
    }

    calc_r_avg();
    V_tot = BP(atoms);

    double r_coeff = 1.0 + (2.0 * uniform() - 1.0) * W_all;

    pos_list atoms_old = atoms;

    for (auto &atom : atoms)
        atom *= r_coeff;

    double dV = BP(atoms) - V_tot;

    if (dV <= 0.0)
    {
        for (auto &atom_spherical : atoms_spherical)
            atom_spherical.x() *= r_coeff;

        r_avg *= r_coeff;
        V_tot += dV;

        return;
    }

    double p = std::min(1.0, std::exp(-beta * dV));

    if (uniform() <= p)
    {
        for (auto &atom_spherical : atoms_spherical)
            atom_spherical.x() *= r_coeff;

        r_avg *= r_coeff;
        V_tot += dV;

        return;
    }

    atoms = std::move(atoms_old);
}

void Fullerene::calc_r_avg()
{
    double sum = 0.0;

    for (auto atom_spherical : atoms_spherical)
        sum += atom_spherical.x();

    r_avg = sum / atoms_spherical.size();
}