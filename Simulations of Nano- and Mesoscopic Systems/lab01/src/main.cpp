#include <iostream>

// Conversion factors
constexpr double Hartree_CF = 27211.386245988;
constexpr double Bohr_CF = 0.0529177210903;

inline double Ha_to_meV(double E) {
    return E * Hartree_CF;
}

inline double meV_to_Ha(double E) {
    return E / Hartree_CF;
}

inline double ab_to_nm(double L) {
    return L * Bohr_CF;
}

inline double nm_to_ab(double L) {
    return L / Bohr_CF;
}

// Functions

class Grid {
public:  
    Grid(size_t n,
            const vec2d_t &delta,
            const vec2d_t &omega,
            double m_eff) : _n(n),
                            _delta(delta),
                            _omega(omega),
                            _m_eff(m_eff),
                            _N_nodes(n * n),
                            _L({delta.first * (n / 2), delta.second * (n / 2)}),
                            _alpha({1.0 / (_m_eff * _omega.first), 1.0 / (_m_eff * _omega.second)}),
                            _node_pos(Eigen::VectorX<vec2d_t>(n * n)),
                            _S{Eigen::MatrixXd::Zero(n * n, n * n)},
                            _K{Eigen::MatrixXd::Zero(n * n, n * n)},
                            _V{Eigen::MatrixXd::Zero(n * n, n * n)},
                            _H{Eigen::MatrixXd::Zero(n * n, n * n)}
    {
        update_node_pos();
        update_matrices();
    }        

private:
    size_t _n;
    vec2d_t _delta;
    vec2d_t _omega;
    double _m_eff;
    size_t _N_nodes;
    vec2d_t _L;
    vec2d_t _alpha;
    Eigen::VectorX<vec2d_t> _node_pos;
    Eigen::MatrixXd _S, _K, _V, _H;

    // void update_node_pos();
    // void update_matrices();
    // double calc_phi(const vec2d_t &pos, const size_t center_node_id) const;
};


// Main function

int main() {

    size_t n = 9;
    double length_x = nm_to_ab(2.0);
    double length_y = nm_to_ab(2.0);
    double omega_x  = meV_to_Ha(80.0);
    double omega_y = meV_to_Ha(200.0);
    double mass_eff = 0.24;

    Grid grid(n, delta, omega, m_eff);

    // zad1
    std::ofstream file_zad1;
    Eigen::Vector3i k_vec{0, 23, 56};

    for (size_t i = 0; i < k_vec.size(); ++i)
    {
        file_zad1.open("data/zad1_" + std::to_string(i) + ".dat");
        grid.print_base_function(k_vec[i], file_zad1, {501, 501});
        file_zad1.close();
    }

    return 0;
}