#include "functions.hpp"

double L = 10;
double x_min = -5.0;
double x_max =  x_min + L;
int nx = 10;
double dx = (x_max - x_min) / (nx - 1);
double y_min = -5.0;
double y_max =  y_min + L;
int ny = 10;
double dy = (y_max - y_min) / (ny - 1);

static auto get_elements_number = [](int new_nx=nx, int new_ny=ny) {return 2*(new_nx-1)*(new_ny-1);};

std::vector<Node> nodes(nx*ny);
std::vector<Element> elements(get_elements_number());

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

double rho(double x, double y) {
    return std::exp( -0.5*(x*x + y*y) );
};

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

void make_grid(int new_nx, int new_ny, double new_x_min, double new_x_max, 
                                       double new_y_min, double new_y_max) {
    if (new_nx != nx || new_ny != ny) {
        nodes.assign(new_nx*new_ny, {});
        elements.assign(get_elements_number(new_nx, new_ny), {});
    }
    nx = new_nx; ny = new_ny;
    x_min = new_x_min; x_max = new_x_max;                                        
    y_min = new_y_min; y_max = new_y_max;                                        
    dx = (x_max - x_min) / (nx - 1);
    dy = (y_max - y_min) / (ny - 1);
    // indeksacja wezlow i elementow zaczyna sie od 0                                   
    for (int j=0; j<ny; ++j)
        for (int i=0; i<nx; ++i) {
            int id = i + j*nx;
            nodes[id] = {id, x_min + i*dx, y_min + j*dy};
        }
    int elem_idx=0;
    for (int j=0; j<ny-1; ++j)
        for (int i=0; i<nx-1; ++i) {
            int id = i + j*nx;
            int n0 = j*nx + i;
            int n1 = n0 + 1;
            int n2 = n0 + nx;
            int n3 = n2 + 1;
            elements[id] = {n0, n1, n2, 0, 1, 2, elem_idx++}; // trojkat 1
            elements[id] = {n1, n2, n3, 0, 1, 2, elem_idx++}; // trojkat 2
        }
}

void shape_functions(double dzeta, double eta, double &phi0, double &phi1, double &phi2) {
    phi0 = -0.5 * (eta + dzeta);
    phi1 =  0.5 * (1 + dzeta);
    phi2 =  0.5 * (1 + eta);
}
// pochodne funkcji ksztaltu po dzeta i eta
static double dphi_dz[3] = { -0.5, 0.5, 0.0 };
static double dphi_de[3] = { -0.5, 0.0, 0.5 };

void jacobian(const Node &A, const Node &B, const Node &C, double dzeta, double eta, double (&J)[2][2], double (&invJ)[2][2], double &detJ) {
    // dzeta, eta nie sa wykorzystywane mimo zaleznosci
    double dx_dz = dphi_dz[0]*A.x + dphi_dz[1]*B.x + dphi_dz[2]*C.x;
    double dx_de = dphi_de[0]*A.x + dphi_de[1]*B.x + dphi_de[2]*C.x;
    double dy_dz = dphi_dz[0]*A.y + dphi_dz[1]*B.y + dphi_dz[2]*C.y;
    double dy_de = dphi_de[0]*A.y + dphi_de[1]*B.y + dphi_de[2]*C.y;

    J[0][0] = dx_dz; J[0][1] = dx_de;
    J[1][0] = dy_dz; J[1][1] = dy_de;

    detJ = J[0][0]*J[1][1] - J[0][1]*J[1][0];

    invJ[0][0] =  J[1][1] / detJ;
    invJ[0][1] = -J[0][1] / detJ;
    invJ[1][0] = -J[1][0] / detJ;
    invJ[1][1] =  J[0][0] / detJ;
}

void assemble_local_matrice_triangle(const std::vector<Node> &nodes, const Element &el, double (&local_matrix)[3][3]) {
    for (int i=0; i<3; ++i) 
        for (int j=0; j<3; ++j) 
            local_matrix[i][j] = 0;

    Node A = nodes[el.i];
    Node B = nodes[el.j];
    Node C = nodes[el.k];

    for (int k=0; k<7; ++k) {
        double dzeta = gauss_points[k][0]; 
        double eta = gauss_points[k][1];
        double w = gauss_weights[k];

        double J[2][2];        // macierz jakobianu
        double invJ[2][2];     // odwrotnosc macierzy jakobianu
        double detJ;           // wyznacznik macierzy jakobianu
        jacobian(A, B, C, dzeta, eta, J, invJ, detJ);

        double dphi_dx[3], dphi_dy[3];
        for (int i=0; i<3; ++i) {
            dphi_dx[i] = invJ[0][0]*dphi_dz[i] + invJ[1][0]*dphi_de[i];
            dphi_dy[i] = invJ[0][1]*dphi_dz[i] + invJ[1][1]*dphi_de[i];
        }

        for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
                local_matrix[i][j] += w * detJ * (dphi_dx[i]*dphi_dx[j] + dphi_dy[i]*dphi_dy[j]);
    }
}

void assemble_local_vector_triangle(const std::vector<Node> &nodes, const Element &el, double (&local_vector)[3]) {
    for (int i=0; i<3; ++i) 
        local_vector[i] = 0.0;

    Node A = nodes[el.i];
    Node B = nodes[el.j];
    Node C = nodes[el.k];

    for (int k=0; k<7; ++k) {
        double dzeta = gauss_points[k][0]; 
        double eta = gauss_points[k][1];
        double w = gauss_weights[k];

        double phi0, phi1, phi2;
        shape_functions(dzeta, eta, phi0, phi1, phi2);
        double phi[3] = {phi0, phi1, phi2};

        double J[2][2];        // macierz jakobianu
        double invJ[2][2];     // odwrotnosc macierzy jakobianu
        double detJ;           // wyznacznik macierzy jakobianu
        jacobian(A, B, C, dzeta, eta, J, invJ, detJ);

        double x = phi0*A.x + phi1*B.x + phi2*C.x;
        double y = phi0*A.y + phi1*B.y + phi2*C.y;
        double rh = rho(x,y);

        for (int i=0; i<3; ++i)
            local_vector[i] += w * detJ * phi[i] * rh;
    }
}

int nr(int local_index, int element_index) {
    int global_index = -1;
    for (auto elem: elements) {
        if (elem.index == element_index) {
            if (elem.a == local_index)
                global_index = elem.i;
            if (elem.b == local_index)
                global_index = elem.j;
            if (elem.c == local_index)
                global_index = elem.k;
        }
    }
    return global_index;
}

// struct Element { 
//     int i, j, k; // globalne indeksy wezlow elementu
//     int a, b, c; // lokalne  indeksy wezlow elementu
//     int index;   // indeks elementu
// };


// bool isBoundaryNode(const Node &p) {
//     return (std::abs(p.x + 5) < 1e-12 ||
//             std::abs(p.x - 5) < 1e-12 ||
//             std::abs(p.y + 5) < 1e-12 ||
//             std::abs(p.y - 5) < 1e-12);
// }