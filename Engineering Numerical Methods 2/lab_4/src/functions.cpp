#include "functions.hpp"

int nx = 3;
double dx = x_max / (nx - 1);
int ny = 3;
double dy = y_max / (ny - 1);

int M = (ny - 1) * (nx - 1);
int num_nodes = nx*ny;
int N_max = 4*nx*ny;
double delta_x = 0.001;

std::vector<Node> global_nodes(num_nodes);
std::vector<std::array<const Node*, 4>> local_nodes(M);
std::vector<std::vector<double>> S(N_max, std::vector<double>(N_max, 0.0));
std::vector<double> F(N_max, 0.0);
std::vector<double> c(N_max, 0.0);

std::vector<std::vector<double>> S_original(N_max, std::vector<double>(N_max, 0.0));
std::vector<double> F_original(N_max, 0.0);

double rho(double x, double y) {
    return std::sin(2*y) * std::sin(x)*std::sin(x);
};

double u_analytical(double x, double y) {
    double exp_2x = std::exp(2*x);
    double exp_2pi = std::exp(2*M_PI);
    return std::sin(2*y) * ( ( (exp_2x*(1.0/exp_2pi - 1) - (exp_2pi - 1))/exp_2x ) / (exp_2pi - 1.0/exp_2pi) + 2.0 - std::cos(2*x)) / 16.0;
};

double diff_u1(const std::function<double(double)>& u, double x, double dx) {
    return (u(x + dx) - u(x - dx)) / (2*dx);
}

double diff_u2(const std::function<double(double)>& u, double x, double dx) {
    return (u(x + dx) - 2*u(x) + u(x - dx)) / (dx*dx);
}

void solve_Ac_b(std::vector<std::vector<double>>& A, std::vector<double>& c, std::vector<double>& b, int N) {
    gsl_matrix *gsl_A = gsl_matrix_alloc(N, N);
    gsl_vector *gsl_b = gsl_vector_alloc(N);
    gsl_vector *gsl_x = gsl_vector_alloc(N);

    for (int i = 0; i < N; ++i) {
        gsl_vector_set(gsl_b, i, b[i]);
        for (int j = 0; j < N; ++j)
            gsl_matrix_set(gsl_A, i, j, A[i][j]);
    }
    // Rozwiązywanie (LU)
    gsl_permutation *p = gsl_permutation_alloc(N);
    int signum;
    gsl_linalg_LU_decomp(gsl_A, p, &signum);
    gsl_linalg_LU_solve(gsl_A, p, gsl_b, gsl_x);

    for (int i = 0; i < N; ++i)
        c[i] = gsl_vector_get(gsl_x, i);

    gsl_matrix_free(gsl_A);
    gsl_vector_free(gsl_b);
    gsl_vector_free(gsl_x);
    gsl_permutation_free(p);
}

double integrate(std::function<double(double)> f, double a, double b, unsigned int n){
    double integral = 0.0;
    double xi, wi;
    gsl_integration_glfixed_table *tab = gsl_integration_glfixed_table_alloc(n);
    for (unsigned int j=0; j<n; ++j) {
        gsl_integration_glfixed_point(a, b, j, &xi, &wi, tab);
        integral += wi*f(xi);
    }
    gsl_integration_glfixed_table_free(tab);
    return integral;
}

auto shape_hermit_functions(int alpha, int i) -> std::function<double(double)> {
    switch(alpha) {
        case 0:
            switch(i) {
                case 0: return [](double xi) { return 0.5 - 0.75*xi + 0.25*xi*xi*xi; };           // funkcja 
                case 1: return [](double xi) { return (1.0 - xi - xi*xi + xi*xi*xi) / 4.0; };     // pochodna
            } break;
        case 1:
            switch(i) {
                case 0: return [](double xi) { return 0.5 + 0.75*xi - 0.25*xi*xi*xi; };           // funkcja
                case 1: return [](double xi) { return (-1.0 - xi + xi*xi + xi*xi*xi) / 4.0; };    // pochodna
            } break;
    }
    return [](double) { return 0.0; };
}

double weight_functions(int i, double xi_1, double xi_2) {
    switch(i) {
        case 0: return 0.25 * (1.0 - xi_1) * (1.0 - xi_2);
        case 1: return 0.25 * (1.0 + xi_1) * (1.0 - xi_2);
        case 2: return 0.25 * (1.0 + xi_1) * (1.0 + xi_2);
        case 3: return 0.25 * (1.0 - xi_1) * (1.0 + xi_2);
    }
    return 0.0;
}

int global_index(int local_index, int phi_i, int phi_j) {
    return 4*(local_index - 1) + phi_i + phi_j*2 + 1;
}

void make_nodes() {
    for (int j=0; j<ny; ++j) {
        for (int i=0; i<nx; ++i) {
            int idx = i + j*nx + 1;
            double x_global = x_min + i*dx;
            double y_global = y_min + j*dy;
            global_nodes[idx-1] = Node{static_cast<unsigned int>(idx), x_global, y_global};

            if (i > 0 && j > 0) {
                int m = i + (j - 1)*(nx - 1);
                int idx1 = i + (j - 1)*nx;
                int idx2 = (i + 1) + (j - 1)*nx;
                int idx3 = (i + 1) + j*nx;
                int idx4 = i + j*nx;

                local_nodes[m-1] = {
                    &global_nodes[idx1-1],
                    &global_nodes[idx2-1],
                    &global_nodes[idx3-1],
                    &global_nodes[idx4-1]
                };
            }
        }
    }
}

double x(int m, double xi_1, double xi_2) {
    double xx = 0.0;
    for (int i=0; i<4; ++i)
        xx += local_nodes[m-1][i]->x * weight_functions(i, xi_1, xi_2);
    return xx;
}

double y(int m, double xi_1, double xi_2) {
    double yy = 0.0;
    for (int i=0; i<4; ++i)
        yy += local_nodes[m-1][i]->y * weight_functions(i, xi_1, xi_2);
    return yy;
}

void update_params(int new_nx, int new_ny) {
    nx = new_nx;
    dx = x_max / (nx - 1);
    ny = new_ny;
    dy = y_max / (ny - 1);

    M = (ny - 1) * (nx - 1);
    num_nodes = nx*ny;
    N_max = 4*nx*ny;

    global_nodes.assign(num_nodes, Node{});
    local_nodes.assign(M, {nullptr, nullptr, nullptr, nullptr});
    S.assign(N_max, std::vector<double>(N_max, 0.0));
    F.assign(N_max, 0.0);
    c.assign(N_max, 0.0);
}

void assemble_global_matrices_MES_2D() {
    for (int m=1; m<=M; ++m) {
        double Jm = (local_nodes[m-1][1]->x - local_nodes[m-1][0]->x) * (local_nodes[m-1][3]->y - local_nodes[m-1][0]->y) / 4.0;
        
        for (int l1=1; l1<=4; ++l1) {
            int alpha1 = l_alpha_beta[l1-1][0];
            int beta1  = l_alpha_beta[l1-1][1];
            for (int i1=0; i1<=1; ++i1) {
                for (int i2=0; i2<=1; ++i2) {
                    // Obliczanie elementu wektora F
                    int p = global_index(local_nodes[m-1][l1-1]->idx, i1, i2);
                    double F_result = integrate(
                        [m, Jm, alpha1, beta1, i1, i2](double xi1) -> double {
                            return integrate(
                                [m, Jm, alpha1, beta1, i1, i2, xi1](double xi2) -> double {
                                    return -Jm * shape_hermit_functions(alpha1, i1)(xi1) * shape_hermit_functions(beta1, i2)(xi2) * rho(x(m, xi1, xi2), y(m, xi1, xi2));
                                }, -1.0, 1.0, 20
                            );
                        }, -1.0, 1.0, 20
                    );
                    F[p-1] += F_result;

                    for (int l2=1; l2<=4; ++l2) {
                        int alpha2 = l_alpha_beta[l2-1][0];
                        int beta2  = l_alpha_beta[l2-1][1];
                        for (int j1=0; j1<=1; ++j1) {
                            for (int j2=0; j2<=1; ++j2) {
                                // Obliczanie elementu macierzy S
                                int q = global_index(local_nodes[m-1][l2-1]->idx, j1, j2);
                                double S_result = integrate(
                                    [Jm, alpha1, beta1, alpha2, beta2, i1, i2, j1, j2](double xi1) -> double {
                                        return integrate(
                                            [Jm, alpha1, beta1, alpha2, beta2, i1, i2, j1, j2, xi1](double xi2) -> double {
                                                return shape_hermit_functions(alpha1, i1)(xi1) * shape_hermit_functions(beta1, i2)(xi2) 
                                                * (diff_u2(shape_hermit_functions(alpha2, j1), xi1) * shape_hermit_functions(beta2, j2)(xi2) 
                                                   + shape_hermit_functions(alpha2, j1)(xi1) * diff_u2(shape_hermit_functions(beta2, j2), xi2)); 
                                            }, -1.0, 1.0, 20
                                        );
                                    }, -1.0, 1.0, 20
                                );
                                S[p-1][q-1] += S_result;
                            }
                        }
                    }

                }
            }
        }
    }

    S_original = S;
    F_original = F;
}

void border_conditions() {
    for (int m=1; m<=M; ++m) {
        for (int l=1; l<=4; ++l) {
            const Node *node = local_nodes[m-1][l-1];
            double eps = 1e-6;
            if (std::abs(node->x - x_min) < eps || std::abs(node->x - x_max) < eps || std::abs(node->y - y_min) < eps || std::abs(node->y - y_max) < eps) {
                for (int i1 = 0; i1 <= 1; ++i1) {
                    for (int i2 = 0; i2 <= 1; ++i2) {
                        if (i1 * i2 == 1)
                            continue;
                        int p = global_index(node->idx, i1, i2);
                        for (int q=0; q<N_max; ++q) {
                            S[p-1][q] = 0.0;
                            S[q][p-1] = 0.0;
                        }
                        S[p-1][p-1] = 1.0;
                        F[p-1] = 0.0;
                    }
                }
            }
        }
    }
}

double functional_integral() {
    double a_num=0.0;
    for (int i=0; i<N_max; ++i) {
        for (int j=0; j<N_max; ++j) {
            a_num -= c[i]*c[j]*S_original[i][j]/2.0; 
        }
        a_num += c[i]*F_original[i];
    }
    return a_num;
}

double u_solution(double x, double y) {
    // Znajdź element, w którym znajduje się punkt (x, y)
    for (int m = 1; m <= M; ++m) {
        double x_m_1 = local_nodes[m-1][0]->x;
        double x_m_2 = local_nodes[m-1][1]->x;
        double y_m_1 = local_nodes[m-1][0]->y;
        double y_m_4 = local_nodes[m-1][3]->y;

        if (x < x_m_1 || x >= x_m_2)
            continue;
        if (y < y_m_1 || y >= y_m_4)
            continue;

        double u = 0.0;
        double xi1 = (x - (x_m_1 + x_m_2) / 2.0) * 2.0/(x_m_2 - x_m_1);
        double xi2 = (y - (y_m_1 + y_m_4) / 2.0) * 2.0/(y_m_4 - y_m_1);

        for (int l = 1; l <= 4; ++l) {
            int alpha = l_alpha_beta[l-1][0];
            int beta  = l_alpha_beta[l-1][1];
            for (int i1 = 0; i1 <= 1; ++i1) {
                for (int i2 = 0; i2 <= 1; ++i2) {
                    int p = global_index(local_nodes[m - 1][l - 1]->idx, i1, i2);
                    u += c[p - 1] * shape_hermit_functions(alpha, i1)(xi1) * shape_hermit_functions(beta, i2)(xi2);
                }
            }
        }
        return u;
    }
    return 0.0;
}