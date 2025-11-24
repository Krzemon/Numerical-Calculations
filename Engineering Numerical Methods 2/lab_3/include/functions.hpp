#pragma once
#include <vector>
#include <cmath>
#include <utility>

inline constexpr double DXI = 0.001;

inline double phi(int i, double xi) {
    switch(i) {
        case 0: return xi*(xi-1)/2.0;
        case 1: return 1.0 - xi*xi;
        case 2: return xi*(xi+1)/2.0;
    }
    return 0.0;
}

/**
 * @brief Oblicza pochodna funkcji bazowej i w punkcie xi metoda roznic centralnych
 */
inline double dphi_num(int i, double xi) {
    double h = DXI;
    return (phi(i, xi+h) - phi(i, xi-h))/(2.0*h);
}

/**
 * @brief Przeksztalca wspolrzedna xi na x w elemencie o koncach xa, xb
 */
inline double x_of_xi(double xa, double xb, double xi) {
    return 0.5*(xa+xb) + 0.5*(xb-xa)*xi;
}

/**
 * @brief Przeksztalca wspolrzedna x na xi w elemencie o koncach xa, xb
 */
inline double xi_of_x(double xa, double xb, double x) {
    return (2.0*x - xa - xb)/(xb - xa);
}

/**
 * @brief Zwraca wezly i wagi kwadratury Gaussa dla 6 punktow
 */
inline void gauss_nodes_weights(std::vector<double>& nodes, std::vector<double>& weights) {
    nodes = {-0.9324695142, -0.6612093864, -0.2386191861,
              0.2386191861, 0.6612093864, 0.9324695142};
    weights = {0.1713244924, 0.3607615730, 0.4679139346,
               0.4679139346, 0.3607615730, 0.1713244924};
}

/**
 * @brief Generuje wezly siatki dla danego M i parametru alpha
 */
inline auto generate_nodes(int M, double alpha, double xmax=6.0) -> std::vector<double> {
    int N = 2*M+1;
    std::vector<double> x(N);
    for(int k=1; k<=N; ++k){
        double arg = (2.0*k - N - 1.0)/N;
        double s = arg>=0 ? 1.0 : -1.0;
        x[k-1] = xmax * std::pow(std::abs(arg), alpha) * s;
    }
    return x;
}

/**
 * @brief Zwraca globalny indeks w macierzach dla elementu m i lokalnego i
 */
inline int global_index(int m, int i) { return 2*m + i; }

/**
 * @brief Składa globalne macierze sztywności S i całek przekrywania O
 */
inline void assemble_global_matrices(int M, const std::vector<double>& xnodes,
                                     std::vector<std::vector<double>>& S,
                                     std::vector<std::vector<double>>& O) {
    int N = 2*M + 1;
    S.assign(N, std::vector<double>(N,0.0));
    O.assign(N, std::vector<double>(N,0.0));

    std::vector<double> g_nodes, g_weights;
    gauss_nodes_weights(g_nodes, g_weights);

    for(int m=0; m<M; ++m){
        double xa = xnodes[2*m];
        double xb = xnodes[2*m+2];
        double Jm = (xb - xa)/2.0;
        for(int i=0;i<3;i++){
            for(int j=0;j<3;j++){
                double Sloc=0.0, Oloc=0.0;
                for(size_t node=0;node<g_nodes.size();node++){
                    double xi = g_nodes[node];
                    double w = g_weights[node];
                    double ph_i = phi(i, xi);
                    double ph_j = phi(j, xi);
                    double dph_i = dphi_num(i, xi);
                    double dph_j = dphi_num(j, xi);
                    double x = x_of_xi(xa, xb, xi);
                    double integr_S = (1.0/(2.0*Jm))*dph_i*dph_j + 0.5*Jm*x*x*ph_i*ph_j;
                    double integr_O = Jm * ph_i*ph_j;
                    Sloc += w*integr_S;
                    Oloc += w*integr_O;
                }
                int p = global_index(m,i);
                int q = global_index(m,j);
                S[p][q] += Sloc;
                O[p][q] += Oloc;
            }
        }
    }

    for(int i=0;i<N;i++)
        for(int j=i+1;j<N;j++){
            double ss = 0.5*(S[i][j]+S[j][i]);
            double oo = 0.5*(O[i][j]+O[j][i]);
            S[i][j]=S[j][i]=ss;
            O[i][j]=O[j][i]=oo;
        }
}