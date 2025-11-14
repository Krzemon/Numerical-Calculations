#pragma once
#include <vector>
#include <cmath>
#include <cassert>
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

inline double dphi_num(int i, double xi) {
    double h = DXI;
    return (phi(i, xi+h) - phi(i, xi-h))/(2.0*h);
}

inline double x_of_xi(double xa, double xb, double xi) {
    return 0.5*(xa+xb) + 0.5*(xb-xa)*xi;
}

inline double xi_of_x(double xa, double xb, double x) {
    return (2.0*x - xa - xb)/(xb - xa);
}

inline void gauss_nodes_weights(std::vector<double>& nodes, std::vector<double>& weights) {
    nodes = {-0.9324695142, -0.6612093864, -0.2386191861,
              0.2386191861, 0.6612093864, 0.9324695142};
    weights = {0.1713244924, 0.3607615730, 0.4679139346,
               0.4679139346, 0.3607615730, 0.1713244924};
}

inline std::vector<double> generate_nodes(int M, double alpha, double xmax=6.0) {
    int N = 2*M+1;
    std::vector<double> x(N);
    for(int k=1; k<=N; ++k){
        double arg = (2.0*k - N - 1.0)/N;
        double s = arg>=0 ? 1.0 : -1.0;
        x[k-1] = xmax * std::pow(std::abs(arg), alpha) * s;
    }
    return x;
}

inline int global_index(int m, int i) { return 2*m + i; }

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
                for(size_t q=0;q<g_nodes.size();q++){
                    double xi = g_nodes[q];
                    double w = g_weights[q];
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
                int qg = global_index(m,j);
                S[p][qg] += Sloc;
                O[p][qg] += Oloc;
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

inline std::vector<std::pair<double,double>> reconstruct_u_on_fine_grid(
    int M, const std::vector<double>& xnodes, const std::vector<double>& c, double dx=0.01)
{
    std::vector<std::pair<double,double>> out;
    double xa = xnodes.front();
    double xb = xnodes.back();

    for(double x=xa; x<=xb+1e-12; x+=dx){
        int m_found=-1;
        for(int m=0;m<M;m++){
            double xa_m=xnodes[2*m], xb_m=xnodes[2*m+2];
            if((x>=xa_m && x<xb_m)||(m==M-1 && x<=xb_m+1e-12)){ m_found=m; break;}
        }
        if(m_found<0) continue;
        double xa_m=xnodes[2*m_found], xb_m=xnodes[2*m_found+2];
        double xi = xi_of_x(xa_m, xb_m, x);
        double u=0.0;
        for(int i=0;i<3;i++){
            int k = global_index(m_found,i);
            u += c[k]*phi(i,xi);
        }
        out.emplace_back(x,u);
    }
    return out;
}