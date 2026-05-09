#include "QPC.hpp"
#include "config.hpp"
#include "params.hpp"

#include <fstream>
#include <cmath>
#include <vector>

QPCSolver::QPCSolver(Grid grid) : m_grid(grid) {}

void QPCSolver::dispersion(const std::string& filename, int N_kx) const
{
    auto file = prepareDataFile(filename);
    stream_config(file.out);

    double a = m_grid.alpha();
    double dy = m_grid.dy();
    double kx_max = M_PI / dy;

    Eigen::VectorXd diag(m_grid.Ny);
    Eigen::VectorXd sub(m_grid.Ny - 1);
    sub.fill(-a);

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver;

    std::vector<double> kx_vals(N_kx);
    std::vector<Eigen::VectorXd> E_vals(N_kx);

    #pragma omp parallel for
    for (int i = 0; i < N_kx; ++i) {
        double kx = -kx_max + (2.0 * kx_max / (N_kx - 1)) * i;
        Eigen::VectorXd diag_local(m_grid.Ny);
        diag_local.fill(4.0 * a - 2.0 * a * std::cos(kx * dy));
        
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> local_solver;
        local_solver.computeFromTridiagonal(diag_local, sub, Eigen::EigenvaluesOnly);
        
        kx_vals[i] = kx;
        E_vals[i] = local_solver.eigenvalues();
    }

    for (int i = 0; i < N_kx; ++i) {
        file.out << kx_vals[i] / bohr_ratio_nm << " " << E_vals[i].transpose() * E_h << "\n";
    }
    
    std::cout << " [saved]: " << file.path << "\n";
}

Modes QPCSolver::modes(double E,
                       const std::string& filename,
                       const std::string& filename_u) const
{
    double dy = m_grid.dy();

    Eigen::MatrixXd tau = make_tau();
    Eigen::MatrixXd H0  = make_H0();

    Eigen::MatrixXd H_gen = Eigen::MatrixXd::Zero(2 * m_grid.Ny, 2 * m_grid.Ny);
    Eigen::MatrixXd S_gen = Eigen::MatrixXd::Zero(2 * m_grid.Ny, 2 * m_grid.Ny);

    H_gen.topRightCorner(m_grid.Ny, m_grid.Ny) = Eigen::MatrixXd::Identity(m_grid.Ny, m_grid.Ny);
    H_gen.bottomLeftCorner(m_grid.Ny, m_grid.Ny) = -tau;
    H_gen.bottomRightCorner(m_grid.Ny, m_grid.Ny) = E * Eigen::MatrixXd::Identity(m_grid.Ny, m_grid.Ny) - H0;

    S_gen.topLeftCorner(m_grid.Ny, m_grid.Ny) = Eigen::MatrixXd::Identity(m_grid.Ny, m_grid.Ny);
    S_gen.bottomRightCorner(m_grid.Ny, m_grid.Ny) = tau.transpose();

    Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> eig_solver;
    eig_solver.compute(H_gen, S_gen);

    Modes result;
    auto fh = filename.empty() ? FileHandle{} : prepareDataFile(filename + ".dat");
    if (!filename.empty()) stream_config(fh.out);

    for (int idx = 0; idx < 2 * m_grid.Ny; ++idx) {
        Eigen::dcomplex lambda = eig_solver.eigenvalues()(idx);
        
        if (std::abs(std::abs(lambda) - 1.0) > 1e-6) continue;

        Eigen::VectorXcd u_mode = eig_solver.eigenvectors().col(idx).head(m_grid.Ny);
        double vel = -2.0 * dy * (lambda * u_mode.dot(tau.transpose() * u_mode)).imag();
        double kx = std::log(lambda).imag() / dy;
        
        // Zapisujemy WSZYSTKIE mody propagujące (zarówno v>0 jak i v<0)
        if (!filename.empty() && std::abs(vel) > 1e-10) {
            fh.out << kx / bohr_ratio_nm << " " << au_to_eV(E) << "\n";
        }

        if (vel > 0.0) {
            result.lambda.push_back(lambda);
            result.u.push_back(u_mode.normalized());
            result.v.push_back(vel);
        }
    }

    if (!filename.empty()) {
        std::cout << " [saved]: " << fh.path << "\n";
    }

    if (!filename_u.empty()) {
        for (size_t n = 0; n < result.u.size(); ++n) {
            auto fh_u = prepareDataFile(filename_u + std::to_string(n) + ".dat");
            stream_config(fh_u.out);
            for (int j = 0; j < m_grid.Ny; ++j) {
                double y = m_grid.y_min + (j + 1) * m_grid.dy();
                fh_u.out << au_to_nm(y) << " "
                         << result.u[n](j).real() << " "
                         << result.u[n](j).imag() << "\n";
            }
            std::cout << " [saved]: " << fh_u.path << "\n";
        }
    }

    return result;
}


std::pair<double,double> QPCSolver::transmission(double E, double V_gates,
                                                  const std::string& filename) const
{
    Modes m = modes(E);
    if (m.lambda.empty()) return {0.0, 0.0};

    Eigen::MatrixXcd M = make_QTBM(m, E, V_gates);
    auto [T, R] = solve_TR(m, M);

    if (!filename.empty()) {
        Eigen::SparseMatrix<Eigen::dcomplex> Msp = M.sparseView();
        Msp.makeCompressed();
        Eigen::SparseLU<Eigen::SparseMatrix<Eigen::dcomplex>> solver;
        solver.analyzePattern(Msp);
        solver.factorize(Msp);

        Eigen::VectorXd density = Eigen::VectorXd::Zero(m_grid.Nx * m_grid.Ny);
        Eigen::VectorXcd b = Eigen::VectorXcd::Zero(m_grid.Nx * m_grid.Ny);

        for (size_t n = 0; n < m.lambda.size(); ++n) {
            Eigen::dcomplex delta_p = 1.0 - 1.0 / m.lambda[n];
            Eigen::dcomplex delta_m = 1.0 - m.lambda[n];
            b.head(m_grid.Ny) = -m_grid.alpha() * m.u[n] * (delta_p - delta_m);

            Eigen::VectorXcd psi = solver.solve(b);
            density += psi.cwiseAbs2();
        }

        density /= static_cast<double>(m.lambda.size());

        auto fh = prepareDataFile(filename + ".dat");
        stream_config(fh.out);
        for (int i = 0; i < m_grid.Nx; ++i) {
            for (int j = 0; j < m_grid.Ny; ++j) {
                double x = m_grid.x_min + (i + 1) * m_grid.dx();
                double y = m_grid.y_min + (j + 1) * m_grid.dy();
                fh.out << au_to_nm(x) << " " << au_to_nm(y) << " "
                       << density(i * m_grid.Ny + j) << "\n";
            }
        }
        std::cout << " [saved]: " << fh.path << "\n";
    }

    return {T, R};
}

void QPCSolver::conductance_vs_gate(double E, double Vg_min, double Vg_max, int steps,
                                    const std::string& filename,
                                    const std::string& filename_V) const
{
    Modes m = modes(E);
    if (m.lambda.empty()) return;

    {
        auto fh = prepareDataFile(filename_V + ".dat");
        stream_config(fh.out);
        double Vg_ref = eV_to_au(-1.0);
        for (int i = 0; i < m_grid.Nx; ++i) {
            for (int j = 0; j < m_grid.Ny; ++j) {
                double x = m_grid.x_min + (i + 1) * m_grid.dx();
                double y = m_grid.y_min + (j + 1) * m_grid.dy();
                fh.out << au_to_nm(x) << " " << au_to_nm(y) << " "
                       << au_to_meV(V_gate(x, y, Vg_ref)) << "\n";
            }
        }
        std::cout << " [saved]: " << fh.path << "\n";
    }

    std::vector<double> Vg_values(steps);
    std::vector<double> T_values(steps);
    std::vector<double> R_values(steps);

    for (int s = 0; s < steps; ++s) {
        Vg_values[s] = Vg_min + s * (Vg_max - Vg_min) / (steps - 1);
    }

    #pragma omp parallel for
    for (int s = 0; s < steps; ++s) {
        double Vg_au = eV_to_au(Vg_values[s]);
        Eigen::MatrixXcd M = make_QTBM(m, E, Vg_au);
        auto [T, R] = solve_TR(m, M);
        T_values[s] = T;
        R_values[s] = R;
    }

    auto fh = prepareDataFile(filename + ".dat");
    stream_config(fh.out);
    fh.out << "# Vgates[eV]  T  R  T+R\n";
    for (int s = 0; s < steps; ++s) {
        fh.out << Vg_values[s] << " " << T_values[s] << " " 
               << R_values[s] << " " << T_values[s] + R_values[s] << "\n";
    }

    std::cout << " [saved]: " << fh.path << "\n";
}

// =====================================================
//  Pomocnicze metody
// =====================================================
Eigen::MatrixXd QPCSolver::make_tau() const
{
    Eigen::MatrixXd tau = Eigen::MatrixXd::Zero(m_grid.Ny, m_grid.Ny);
    tau.diagonal().fill(-m_grid.alpha());
    return tau;
}

Eigen::MatrixXd QPCSolver::make_H0(double E) const
{
    double a = m_grid.alpha();
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(m_grid.Ny, m_grid.Ny);
    H.diagonal(-1).fill(-a);
    H.diagonal(0).fill(4.0 * a);
    H.diagonal(1).fill(-a);
    if (E != 0.0) H.diagonal(0).array() -= E;
    return H;
}

double QPCSolver::V_gate(double x, double y, double V_gates) const
{
    double sigma = nm_to_au(QPCParams::sigma_nm);
    double arg1 = (x * x) / (sigma * sigma) + 
                  (y - m_grid.y_min) * (y - m_grid.y_min) / (sigma * sigma);
    double arg2 = (x * x) / (sigma * sigma) + 
                  (y - m_grid.y_max) * (y - m_grid.y_max) / (sigma * sigma);
    
    return QPCParams::V_prefactor * V_gates * 
           (std::exp(-arg1 * arg1) + std::exp(-arg2 * arg2));
}

Eigen::MatrixXcd QPCSolver::make_QTBM(const Modes& m, double E, double V_gates) const
{
    double dx = m_grid.dx();
    double dy = m_grid.dy();
    int N = m_grid.Nx * m_grid.Ny;

    Eigen::MatrixXd tau = make_tau();
    Eigen::MatrixXd H0  = make_H0(E);

    Eigen::MatrixXcd M = Eigen::MatrixXcd::Zero(N, N);

    for (int k = 0; k < m_grid.Nx; ++k) {
        int offset = k * m_grid.Ny;
        M.block(offset, offset, m_grid.Ny, m_grid.Ny) = H0;

        if (k + 1 < m_grid.Nx) {
            M.block(offset, offset + m_grid.Ny, m_grid.Ny, m_grid.Ny) = tau.transpose();
            M.block(offset + m_grid.Ny, offset, m_grid.Ny, m_grid.Ny) = tau;
        }
    }

    if (V_gates != 0.0) {
        for (int i = 0; i < m_grid.Nx; ++i) {
            for (int j = 0; j < m_grid.Ny; ++j) {
                double x = m_grid.x_min + (i + 1) * dx;
                double y = m_grid.y_min + (j + 1) * dy;
                M(i * m_grid.Ny + j, i * m_grid.Ny + j) += V_gate(x, y, V_gates);
            }
        }
    }

    Eigen::MatrixXcd ALPHA = Eigen::MatrixXcd::Zero(m_grid.Ny, m_grid.Ny);
    Eigen::MatrixXcd BETA  = Eigen::MatrixXcd::Zero(m_grid.Ny, m_grid.Ny);

    for (int i = 0; i < m_grid.Ny; ++i) {
        for (int j = 0; j < m_grid.Ny; ++j) {
            for (size_t n = 0; n < m.lambda.size(); ++n) {
                Eigen::dcomplex coeff = std::conj(m.u[n](j)) * m.u[n](i);
                ALPHA(i, j) += coeff * (1.0 - m.lambda[n]);
                BETA(i, j)  += coeff * (1.0 - m.lambda[n]);
            }
        }
    }

    M.topLeftCorner(m_grid.Ny, m_grid.Ny) += 
        tau - tau * ALPHA;
    M.bottomRightCorner(m_grid.Ny, m_grid.Ny) += 
        tau.transpose() - tau.transpose() * BETA;

    return M;
}

std::pair<double,double> QPCSolver::solve_TR(const Modes& m,
                                              const Eigen::MatrixXcd& M) const
{
    Eigen::SparseMatrix<Eigen::dcomplex> M_sparse = M.sparseView();
    M_sparse.makeCompressed();
    Eigen::SparseLU<Eigen::SparseMatrix<Eigen::dcomplex>> lu_solver;
    lu_solver.analyzePattern(M_sparse);
    lu_solver.factorize(M_sparse);

    double T_total = 0.0, R_total = 0.0;
    Eigen::VectorXcd rhs = Eigen::VectorXcd::Zero(m_grid.Nx * m_grid.Ny);

    for (size_t n_in = 0; n_in < m.lambda.size(); ++n_in) {
        Eigen::dcomplex delta_p = 1.0 - 1.0 / m.lambda[n_in];
        Eigen::dcomplex delta_m = 1.0 - m.lambda[n_in];
        rhs.head(m_grid.Ny) = -m_grid.alpha() * m.u[n_in] * (delta_p - delta_m);

        Eigen::VectorXcd psi = lu_solver.solve(rhs);

        for (size_t n_out = 0; n_out < m.lambda.size(); ++n_out) {
            Eigen::dcomplex c_refl = (n_out == n_in ? -1.0 : 0.0);
            Eigen::dcomplex c_trans = 0.0;

            for (int j = 0; j < m_grid.Ny; ++j) {
                c_refl += std::conj(m.u[n_out](j)) * psi(j);
                c_trans += std::conj(m.u[n_out](j)) * psi.tail(m_grid.Ny)(j);
            }

            T_total += std::norm(c_trans) * std::abs(m.v[n_out] / m.v[n_in]);
            R_total += std::norm(c_refl) * std::abs(m.v[n_out] / m.v[n_in]);
        }
    }

    return {T_total, R_total};
}
