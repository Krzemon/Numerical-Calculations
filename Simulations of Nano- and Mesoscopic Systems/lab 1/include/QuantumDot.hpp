#pragma once
#include "params.hpp"
#include <vector>
#include <cmath>
#include <ostream>
#include <stdexcept>
#include <string>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
 
/**
 * @defgroup zad_1 Zadanie 1
 * @defgroup zad_2 Zadanie 2
 * @defgroup zad_3 Zadanie 3
 * @defgroup zad_4 Zadanie 4
 * @defgroup zad_5 Zadanie 5
 */

//   true  -> nanometry [nm]
//   false -> jednostki atomowe [a.u.]
constexpr bool USE_NM = true;
 
class QuantumDot {
public:
    const int    n_grid;   // wezly na kierunek
    const int    N;        // laczna liczba wezlow N = n_grid^2
    const double ddx;      // odstep miedzy wezlami [a.u.]
    const double a;        // polowa szerokosci siatki [-a, a] [a.u.]
    const double alpha_x;  // szerokosc gaussjana w x [a.u.], alpha_x = hbar/(m*wx)
    const double alpha_y;  // szerokosc gaussjana w y [a.u.], alpha_y = hbar/(m*wy)
    const double m_eff;    // masa efektywna [a.u.]
    const double wx;       // czestosc kolowa w x [a.u.]
    const double wy;       // czestosc kolowa w y [a.u.]
 
    // siatka wezlow [a.u.], mapowanie: k = i*n_grid + j
    std::vector<double> xk, yk;
 
    Eigen::MatrixXd H, S;
    Eigen::VectorXd eigenvalues;
    Eigen::MatrixXd eigenvectors;
 
    QuantumDot(int n_nodes = ::n, double dx_nm = 2.0, double m = m_effect,
               double wx_in = omega_x, double wy_in = omega_y) :
        n_grid  (n_nodes),
        N       (n_nodes * n_nodes),
        ddx     (dx_nm / a_0),
        a       (ddx * (n_nodes - 1) / 2.0),
        alpha_x (h_bar / (m * wx_in)),
        alpha_y (h_bar / (m * wy_in)),
        m_eff   (m),
        wx      (wx_in),
        wy      (wy_in),
        H(n_nodes * n_nodes, n_nodes * n_nodes),
        S(n_nodes * n_nodes, n_nodes * n_nodes)
    {
        xk.resize(N);
        yk.resize(N);
        for (int k = 0; k < N; ++k) {
            xk[k] = -a + ddx * (k / n_nodes);
            yk[k] = -a + ddx * (k % n_nodes);
        }
        buildMatrices();
    }
 
    /**
     * @brief Rozwiazanie uogolnionego problemu wlasnego H c = E S c.
     */
    void solve()
    {
        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(H, S);
        if (es.info() != Eigen::Success)
            throw std::runtime_error("Eigen: nie udalo sie rozwiazac problemu wlasnego");
        eigenvalues  = es.eigenvalues();
        eigenvectors = es.eigenvectors();
    }
 
    /**
     * @brief k-ta funkcja bazowa jako iloczyn gaussjanow 1D:
     *        phi_k(x,y) = g(x, xk, alpha_x) * g(y, yk, alpha_y)
     *        g(t, t0, alpha) = (1/alpha/pi)^(1/4) * exp(-(t-t0)^2 / (2*alpha))
     */
    double phi(int k, double x, double y) const
    {
        auto g = [](double t, double t0, double alpha) {
            return std::exp(-std::pow(t - t0, 2) / (2.0 * alpha))/ std::pow(alpha * M_PI, 0.25);
        };
        return g(x, xk[k], alpha_x) * g(y, yk[k], alpha_y);
    }
 
    /**
     * @brief Energia state-tego stanu w eV.
     * @attention wymaga wczesniejszego solve()
     */
    double energy_eV(int state) const
    {
        checkSolved();
        return eigenvalues(state) * E_h;
    }
 
    /**
     * @brief k-ta funkcja bazowa na siatce nx*ny.
     *        Kolumny: x  y  phi_k(x,y)
     * @addtogroup zad_1
     * @attention nie wymaga solve()
     */
    void write_basis(std::ostream& out, int k, int nx = 101, int ny = 101) const
    {
        out << "# units=" << (USE_NM ? "nm" : "au") << "\n";
        const double dx = 2.0 * a / (nx - 1);
        const double dy = 2.0 * a / (ny - 1);
        for (int ix = 0; ix < nx; ++ix) {
            double x = -a + ix * dx;
            for (int iy = 0; iy < ny; ++iy) {
                double y = -a + iy * dy;
                out << toOut(x) << '\t' << toOut(y) << '\t'
                    << phi(k, x, y) << '\n';
            }
            out << '\n';
        }
    }
 
    /**
     * @brief n_states najnizszych energii.
     *        Kolumny: i  E[eV]
     * @addtogroup zad_2
     * @attention wymaga wczesniejszego solve()
     */
    void write_energies(std::ostream& out, int n_states = 10) const
    {
        checkSolved();
        out << "# i\tE[eV]\n";
        for (int i = 0; i < n_states; ++i)
            out << i << '\t' << energy_eV(i) << '\n';
    }

    /**
     * @brief Zapis macierzy (H lub S) do pliku w formacie wierszowym.
     *        Kolumny: i  j  M(i,j)
     * @addtogroup zad_2
     * @attention nie wymaga solve()
     */
    void write_matrix(std::ostream& out, const Eigen::MatrixXd& M, const std::string& name) const
    {
        out << "# matrix=" << name << "  size=" << N << "x" << N << "\n";
        out << "# i\tj\tM(i,j)\n";
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                out << i << '\t' << j << '\t' << M(i, j) << '\n';
    }

    /**
     * @brief mapa funkcji falowej na siatce nx*ny.
     *        Kolumny: x  y  psi  psi^2
     * @addtogroup zad_3
     * @attention wymaga wczesniejszego solve()
     */
    void write_wavefunction(std::ostream& out, int state, int nx = 101, int ny = 101) const
    {
        checkSolved();
        out << "# units=" << (USE_NM ? "nm" : "au")
            << "  state=" << state << "  E=" << energy_eV(state) << " eV\n";
        const double dx = 2.0 * a / (nx - 1);
        const double dy = 2.0 * a / (ny - 1);
        for (int ix = 0; ix < nx; ++ix) {
            double x = -a + ix * dx;
            for (int iy = 0; iy < ny; ++iy) {
                double y   = -a + iy * dy;
                double psi = 0.0;
                for (int k = 0; k < N; ++k)
                    psi += eigenvectors(k, state) * phi(k, x, y);
                out << toOut(x) << '\t' << toOut(y) << '\t'
                    << psi << '\t' << psi * psi << '\n';
            }
            out << '\n';
        }
    }
 
    /**
     * @brief naglowek tabeli E(wx).
     *        Format: hbar_wx[meV]  E_0[eV]  E_1[eV] ...
     * @addtogroup zad_4
     */
    static void write_energy_vs_omega_header(std::ostream& out, int n_states)
    {
        out << "# hbar_wx[meV]";
        for (int s = 0; s < n_states; ++s) out << "\tE_" << s << "[eV]";
        out << '\n';
    }
 
    /**
     * @brief jeden wiersz tabeli E(wx).
     * @addtogroup zad_4
     * @attention wymaga wczesniejszego solve()
     */
    void write_energy_vs_omega_row(std::ostream& out, double hbwx_meV, int n_states) const
    {
        checkSolved();
        out << hbwx_meV;
        for (int s = 0; s < n_states; ++s) out << '\t' << energy_eV(s);
        out << '\n';
    }
 
private:
    /**
     * @brief Konwersja jednostek z a.u. na nm
     */
    static double toOut(double x_au)
    {
        return USE_NM ? x_au * a_0 : x_au;
    }
 
    /**
     * @brief Obliczenie macierzy H i S dla bazy gaussjanow.
     * @addtogroup zad_2
     */
    void buildMatrices()
    {
        for (int k = 0; k < N; ++k) {
            for (int l = 0; l < N; ++l) {
                const double x_diff = xk[k] - xk[l];
                const double y_diff = yk[k] - yk[l];
                const double x_sum  = xk[k] + xk[l];
                const double y_sum  = yk[k] + yk[l];
 
                const double S_kl = std::exp(
                    - x_diff * x_diff / (4.0 * alpha_x)
                    - y_diff * y_diff / (4.0 * alpha_y)
                );
 
                const double K_kl = -0.5 / m_eff * (
                    (x_diff * x_diff - 2.0 * alpha_x) / (4.0 * alpha_x * alpha_x) +
                    (y_diff * y_diff - 2.0 * alpha_y) / (4.0 * alpha_y * alpha_y)
                ) * S_kl;
 
                const double V_kl = 0.5 * m_eff * (
                    wx * wx * (x_sum * x_sum + 2.0 * alpha_x) / 4.0 +
                    wy * wy * (y_sum * y_sum + 2.0 * alpha_y) / 4.0
                ) * S_kl;
 
                S(k, l) = S_kl;
                H(k, l) = K_kl + V_kl;
            }
        }
    }
 
    void checkSolved() const
    {
        if (eigenvalues.size() == 0)
            throw std::runtime_error(
                "QuantumDot::solve() nie zostalo wywolane przed dostepem do wynikow.");
    }
};
 