#include "functions.hpp"

// ------------------
// PARAMETRY GLOBALNE
// ------------------
int g_nx = 0, g_ny = 0;
double g_xmin = 0, g_xmax = 0, g_ymin = 0, g_ymax = 0;
double g_dx = 0, g_dy = 0;
int g_N = 0, g_M = 0;

std::vector<Node> g_nodes;
std::vector<NodeSet> g_elem_nodes;
std::vector<double> g_E;
std::vector<double> g_O;
std::vector<std::array<double,9>> g_E_local;
std::vector<std::array<double,9>> g_O_local;
std::vector<double> evals;
std::vector<std::vector<double>> evecs;

static const double gauss_points[7][2] = {
    { -0.333333333333333, -0.333333333333333 },
    { -0.059715871789770, -0.059715871789770 },
    { -0.059715871789770, -0.880568256420460 },
    { -0.880568256420460, -0.059715871789770 },
    { -0.797426985353088, -0.797426985353088 },
    { -0.797426985353088,  0.594853970706174 },
    {  0.594853970706174, -0.797426985353088 }
};

static const double gauss_weights[7] = {
    0.450000000000000,
    0.264788305577012,
    0.264788305577012,
    0.264788305577012,
    0.251878361089654,
    0.251878361089654,
    0.251878361089654
};

const int lg(int elem_idx, int local_node_idx) {
    return g_elem_nodes[elem_idx][local_node_idx]->idx;
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

double integrate_over_reference_space(std::function<double(double,double)> f) {
    double sum = 0.0;
    for(int k=0; k<7; ++k) 
        sum += gauss_weights[k] * f(gauss_points[k][0],gauss_points[k][1]);
    return sum;
}

double x_map(int elem_idx, double dzeta, double eta){
    double x=0.0;
    for(int i=0; i<3; ++i) 
        x += g_elem_nodes[elem_idx][i]->x * phi[i](dzeta,eta);
    return x;
}

double y_map(int elem_idx, double dzeta, double eta){
    double y=0.0;
    for(int i=0; i<3; ++i) 
        y += g_elem_nodes[elem_idx][i]->y * phi[i](dzeta,eta);
    return y;
}

auto static xm_fn(int elem_idx) -> std::function<double(double,double)> {
    return [elem_idx](double dzeta, double eta) { return x_map(elem_idx, dzeta, eta); };
}
auto static ym_fn(int elem_idx) -> std::function<double(double,double)> {
    return [elem_idx](double dzeta, double eta) { return y_map(elem_idx, dzeta, eta); };
}

double jacobian(int elem_idx, double dzeta, double eta){
    double dx_ddzeta = df_num_2D_x(xm_fn(elem_idx), dzeta, eta);
    double dx_deta   = df_num_2D_y(xm_fn(elem_idx), dzeta, eta);
    double dy_ddzeta = df_num_2D_x(ym_fn(elem_idx), dzeta, eta);
    double dy_deta   = df_num_2D_y(ym_fn(elem_idx), dzeta, eta);
    return dx_ddzeta * dy_deta - dx_deta * dy_ddzeta;
}

void make_grid(int nx, double x_min, double x_max, int ny, double y_min, double y_max){
    g_nx=nx; g_ny=ny; g_xmin=x_min; g_xmax=x_max; g_ymin=y_min; g_ymax=y_max;
    g_dx=(x_max-x_min)/(nx-1);
    g_dy=(y_max-y_min)/(ny-1);
    g_N=nx*ny;
    g_M=2*(nx-1)*(ny-1);

    g_nodes.resize(g_N);
    g_elem_nodes.resize(g_M);
    g_E_local.assign(g_M,std::array<double,9>{});
    g_O_local.assign(g_M,std::array<double,9>{});

    g_E.assign(g_N*g_N,0.0);
    g_O.assign(g_N*g_N,0.0);

    fill_node_data();
    fill_local_matrices();
    assemble_matrices();
    apply_boundary_conditions();
}

void fill_node_data(){
    for(int i=0; i<g_nx; ++i){
        for(int j=0; j<g_ny; ++j){
            double x = g_xmin + i*g_dx;
            double y = g_ymin + j*g_dy;
            int idx = i + j*g_nx;
            g_nodes[idx] = {x,y,idx};
            if(i>0 && j>0){
                int m = 2*((i-1) + (j-1)*(g_nx-1));
                int idx0 = (i-1)+(j-1)*g_nx;
                int idx1 = i + (j-1)*g_nx;
                int idx2 = i + j*g_nx;
                int idx3 = (i-1) + j*g_nx;

                g_elem_nodes[m] = { &g_nodes[idx0], &g_nodes[idx1], &g_nodes[idx3] };
                g_elem_nodes[m+1] = { &g_nodes[idx2], &g_nodes[idx3], &g_nodes[idx1] };
            }
        }
    }
}

void fill_local_matrices(){
    for(int m=0; m<g_M; ++m){
        for(int i=0; i<3; ++i){
            for(int j=0; j<3; ++j){
                auto E_integrand = [m,i,j](double dzeta, double eta){
                    double J = jacobian(m, dzeta, eta);
                    double J_det = J;

                    // x(dzeta,eta), y(dzeta,eta)
                    auto xm = xm_fn(m);
                    auto ym = ym_fn(m);

                    double dx_ddzeta  = df_num_2D_x(xm, dzeta, eta);
                    double dx_deta = df_num_2D_y(xm, dzeta, eta);
                    double dy_ddzeta  = df_num_2D_x(ym, dzeta, eta);
                    double dy_deta = df_num_2D_y(ym, dzeta, eta);

                    // macierz odwrotna J^(-1)
                    double invJ[2][2];
                    invJ[0][0] =  dy_deta / J_det;
                    invJ[0][1] = -dx_deta / J_det;
                    invJ[1][0] = -dy_ddzeta  / J_det;
                    invJ[1][1] =  dx_ddzeta  / J_det;

                    // pochodne funkcji kształtu phi_i
                    double dphi_i_dxi  = df_num_2D_x(phi[i], dzeta, eta);
                    double dphi_i_deta = df_num_2D_y(phi[i], dzeta, eta);
                    double dphi_j_dxi  = df_num_2D_x(phi[j], dzeta, eta);
                    double dphi_j_deta = df_num_2D_y(phi[j], dzeta, eta);

                    double grad_i_x = dphi_i_dxi*invJ[0][0] + dphi_i_deta*invJ[1][0];
                    double grad_i_y = dphi_i_dxi*invJ[0][1] + dphi_i_deta*invJ[1][1];
                    double grad_j_x = dphi_j_dxi*invJ[0][0] + dphi_j_deta*invJ[1][0];
                    double grad_j_y = dphi_j_dxi*invJ[0][1] + dphi_j_deta*invJ[1][1];

                    return J * (grad_i_x*grad_j_x + grad_i_y*grad_j_y);
                };

                g_E_local[m][i*3 + j] = integrate_over_reference_space(E_integrand);

                auto O_integrand = [m,i,j](double dzeta, double eta){
                    return jacobian(m, dzeta, eta) * phi[i](dzeta, eta) * phi[j](dzeta, eta);
                };

                g_O_local[m][i*3 + j] = integrate_over_reference_space(O_integrand);
            }
        }
    }
}

void assemble_matrices(){
    for(int m=0; m<g_M; ++m){
        for(int i=0; i<3; ++i){
            // int p = g_elem_nodes[m][i]->idx;
            int p = lg(m,i);
            for(int j=0; j<3; ++j){
                // int q = g_elem_nodes[m][j]->idx;
                int q = lg(m,j);
                g_E[p*g_N + q] += g_E_local[m][i*3 + j];
                g_O[p*g_N + q] += g_O_local[m][i*3 + j];
            }
        }
    }
}

void apply_boundary_conditions(){
    double eps = 1e-6;
    for(int i=0 ;i<g_N; ++i){
        const Node &node = g_nodes[i];
        if(std::abs(node.x - g_xmin)<eps || std::abs(node.x - g_xmax)<eps ||
           std::abs(node.y - g_ymin)<eps || std::abs(node.y - g_ymax)<eps) {
            for(int j=0; j<g_N; ++j){
                g_E[i*g_N + j] = 0.0; 
                g_E[j*g_N + i] = 0.0;
                g_O[i*g_N + j] = 0.0; 
                g_O[j*g_N + i] = 0.0;
            }
            g_E[i*g_N + i] = 2000.0;
            g_O[i*g_N + i] = 1.0;
        }
    }
}

// void newmark_time_integration(double dt, int n_steps, const std::vector<double>& u0, std::vector<double>& u_curr){
//     std::vector<double> u_prev = u0;
//     u_curr = u0;
//     std::vector<double> u_next(g_N,0.0);
//     std::vector<double> RHS(g_N,0.0);

//     for(int n=0;n<n_steps;n++){
//         for(int i=0;i<g_N;i++){
//             RHS[i]=0.0;
//             for(int j=0;j<g_N;j++){
//                 RHS[i] += 2.0*g_O[i*g_N+j]*u_curr[j] - g_O[i*g_N+j]*u_prev[j] - dt*dt*g_E[i*g_N+j]*u_curr[j];
//             }
//         }
//         u_next = solve_linear_system(g_O,RHS);
//         u_prev = u_curr;
//         u_curr = u_next;
//     }
// }

std::vector<double> solve_linear_system(const std::vector<double>& A, const std::vector<double>& b){
    int n = b.size();
    gsl_matrix* gA = gsl_matrix_alloc(n,n);
    gsl_vector* gb = gsl_vector_alloc(n);
    gsl_vector* gx = gsl_vector_alloc(n);

    for(int i=0;i<n;i++){
        gsl_vector_set(gb,i,b[i]);
        for(int j=0;j<n;j++) gsl_matrix_set(gA,i,j,A[i*n+j]);
    }

    gsl_permutation* p = gsl_permutation_alloc(n);
    int signum;
    gsl_linalg_LU_decomp(gA,p,&signum);
    gsl_linalg_LU_solve(gA,p,gb,gx);

    std::vector<double> x(n);
    for(int i=0;i<n;i++) x[i]=gsl_vector_get(gx,i);

    gsl_permutation_free(p);
    gsl_matrix_free(gA);
    gsl_vector_free(gb);
    gsl_vector_free(gx);
    return x;
}

/**
 * @brief Rozwiazuje problem wlasny z wykorzystaniem lapack
 */
void solve_generalized_eigen_lapack(const double* E, const double* O, int N,
                                    std::vector<double>& evals,
                                    std::vector<std::vector<double>>& evecs) {
    std::vector<double> A(N*N), B(N*N);

    // Zamiana na kolumnowy układ LAPACK
    for(int i = 0; i < N; ++i)
        for(int j = 0; j < N; ++j){
            A[j*N + i] = E[i*N + j];
            B[j*N + i] = O[i*N + j];
        }

    std::vector<double> w(N);

    int itype = 1;
    char jobz = 'V';
    char uplo = 'U';
    int lda = N, ldb = N;
    int info;

    // --- Query buforów roboczych ---
    int lwork = -1, liwork = -1;
    double work_query;
    int iwork_query;

    dsygvd_(&itype, &jobz, &uplo, &N,
            A.data(), &lda, B.data(), &ldb,
            w.data(), &work_query, &lwork, &iwork_query, &liwork, &info);

    if(info != 0){
        std::cerr << "dsygvd_ (query) returned " << info << "\n";
        return;
    }

    lwork = static_cast<int>(work_query);
    liwork = iwork_query;

    std::vector<double> work(lwork);
    std::vector<int> iwork(liwork);

    // --- Właściwe wywołanie LAPACK ---
    dsygvd_(&itype, &jobz, &uplo, &N,
            A.data(), &lda, B.data(), &ldb,
            w.data(), work.data(), &lwork, iwork.data(), &liwork, &info);

    if(info){
        std::cerr << "dsygvd_ returned " << info << "\n";
        return;
    }

    // Kopiowanie wyników
    evals.assign(N, 0.0);
    evecs.assign(N, std::vector<double>(N, 0.0));

    for(int k = 0; k < N; ++k){
        evals[k] = w[k];
        for(int i = 0; i < N; ++i){
            evecs[i][k] = A[k*N + i]; // kolumny to wektory własne
        }
    }
}