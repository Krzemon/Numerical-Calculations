#include "solution.hpp"
#include "config.hpp"

void calc_zad_1(std::ofstream& out) {
    make_grid(10, -5, 5, 10, -5, 5);
    out << std::fixed << std::setprecision(6);
    out << "# index x y nx=" << g_nx << " ny=" << g_ny << "\n";
    for (const auto& node : g_nodes) {
        out << node.idx << " " << node.x << " " << node.y << "\n";
    }

    out << "\n# elementy i j k\n";
    for (std::size_t m = 0; m < g_elem_nodes.size(); ++m) {
        const NodeSet& e = g_elem_nodes[m];
        out << m << " "
            << e[0]->idx << " "
            << e[1]->idx << " "
            << e[2]->idx << "\n";
    }
}

void save_global_E(std::ofstream& out) {
    out << std::fixed << std::setprecision(6);
    for(int i=0; i<g_N; ++i){
        for(int j=0; j<g_N; ++j){
            out << g_E[i*g_N + j] << " ";
        }
        out << "\n";
    }
}

void save_global_O(std::ofstream& out) {
    out << std::fixed << std::setprecision(6);
    for(int i=0; i<g_N; ++i){
        for(int j=0; j<g_N; ++j){
            out << g_O[i*g_N + j] << " ";
        }
        out << "\n";
    }
}

void calc_zad_2(std::ofstream &E_out, std::ofstream &O_out, std::ofstream &evals_out, std::ofstream &modes_out) {
    save_global_E(E_out);
    save_global_O(O_out);
    solve_generalized_eigen_lapack(g_E.data(), g_O.data(), g_N, evals, evecs);
    
    int K = std::min(10, g_N);

    evals_out << std::setprecision(15);
    for(int k = 0; k < K; k++){
        double lambda = evals[k];
        if(lambda < 0 && lambda > -1e-12) lambda = 0.0;
        evals_out << lambda << "\n";
    }

    const int Nx = 200;
    const int Ny = 200;

    for(int k = 0; k < K; ++k){
        std::string fname = "mode_" + std::to_string(k) + ".dat";
        auto fh = prepareDataFile(fname);  // przygotowanie pliku
        std::ofstream &fout = fh.out;

        fout << std::setprecision(12);
        fout << "# lambda " << evals[k] << "\n";
        fout << "# x y phi\n";

        double dx = (g_xmax - g_xmin) / (Nx - 1);
        double dy = (g_ymax - g_ymin) / (Ny - 1);

        for(int iy = 0; iy < Ny; ++iy){
            double y = g_ymin + iy*dy;
            for(int ix = 0; ix < Nx; ++ix){
                double x = g_xmin + ix*dx;

                double uk = 0.0;
                bool found = false;

                // iterujemy po wszystkich elementach (dwa trojkaty na kwadrat)
                for(int m = 0; m < g_M; m += 2){
                    double x_left   = g_elem_nodes[m][0]->x;
                    double x_right  = g_elem_nodes[m][1]->x;
                    double y_bottom = g_elem_nodes[m][0]->y;
                    double y_top    = g_elem_nodes[m][2]->y;

                    if(x < x_left || x > x_right || y < y_bottom || y > y_top)
                        continue;

                    // wybór trójkąta w kwadracie
                    int elem_idx = ((x - x_left) + (y - y_bottom) > (x_right - x_left)) ? m + 1 : m;
                    const NodeSet &tri = g_elem_nodes[elem_idx];

                    // mały układ liniowy do współrzędnych lokalnych
                    double x0 = tri[0]->x, y0 = tri[0]->y;
                    double x1 = tri[1]->x, y1 = tri[1]->y;
                    double x2 = tri[2]->x, y2 = tri[2]->y;

                    double A[2][2] = { { x1 - x0, x2 - x0 }, { y1 - y0, y2 - y0 } };
                    double b[2] = { x - x0, y - y0 };

                    // rozwiąż 2x2 układ: A * xi = b
                    double det = A[0][0]*A[1][1] - A[0][1]*A[1][0];
                    if(std::fabs(det) < 1e-14) continue;
                    double xi1 = ( b[0]*A[1][1] - b[1]*A[0][1] ) / det;
                    double xi2 = ( A[0][0]*b[1] - A[1][0]*b[0] ) / det;

                    // funkcje kształtu dla trójkąta liniowego
                    double phi0 = 1.0 - xi1 - xi2;
                    double phi1 = xi1;
                    double phi2 = xi2;

                    uk = phi0 * evecs[tri[0]->idx][k] +
                        phi1 * evecs[tri[1]->idx][k] +
                        phi2 * evecs[tri[2]->idx][k];

                    found = true;
                    break;
                }

                if(!found) uk = 0.0;  // punkt poza elementami
                fout << x << " " << y << " " << uk << "\n";
            }
        }
        std::cout << "Zapisano: " << fh.path << "\n";
    }

}


void calc_zad_3(std::ofstream& out){
    // std::vector<double> eigvals;
    // std::vector<std::vector<double>> modes;
    // solve_eigenproblem(eigvals,modes);

    // if(modes.size() < 2) throw std::runtime_error("Za mało stanów własnych!");
    // auto ccc2 = modes[1]; // drugi mod
    // out << std::fixed << std::setprecision(6);
    // for(double val : ccc2) out << val << "\n";
}

// ------------------
// ZADANIE 4: zapis przykładowego stanu własnego
// ------------------
void calc_zad_4(std::ofstream& out){
    // std::vector<double> eigvals;
    // std::vector<std::vector<double>> modes;
    // solve_eigenproblem(eigvals,modes);

    // if(modes.size() < 2) throw std::runtime_error("Za mało stanów własnych!");
    // auto y0 = modes[1];

    // out << std::fixed << std::setprecision(6);
    // for(double val : y0) out << val << "\n";
}

// ------------------
// ZADANIE 5: integracja w czasie metodą Newmarka
// ------------------
void calc_zad_5(std::ofstream& out){
    // std::vector<double> eigvals;
    // std::vector<std::vector<double>> modes;
    // solve_eigenproblem(eigvals,modes);

    // if(modes.size() < 2) throw std::runtime_error("Za mało stanów własnych!");
    // auto u0 = modes[1]; // startowy wektor
    // std::vector<double> u_curr;
    // double dt = 0.01;
    // int n_steps = 50;

    // newmark_time_integration(dt,n_steps,u0,u_curr);

    // out << std::fixed << std::setprecision(6);
    // for(double val : u_curr) out << val << "\n";
}

// ------------------
// ZADANIE 6: mapowanie rozwiązania u(x,y,t) i zapis map
// ------------------
void calc_zad_6(std::ofstream& out){
    // std::vector<double> eigvals;
    // std::vector<std::vector<double>> modes;
    // solve_eigenproblem(eigvals,modes);

    // if(modes.size() < 2) throw std::runtime_error("Za mało stanów własnych!");
    // auto yyyk = modes[1]; // przykładowy wektor współczynników

    // int nx_map = 20;
    // int ny_map = 20;
    // double dx_map = (g_xmax - g_xmin)/(nx_map-1);
    // double dy_map = (g_ymax - g_ymin)/(ny_map-1);

    // out << std::fixed << std::setprecision(6);
    // for(int j=0;j<ny_map;j++){
    //     double y = g_ymin + j*dy_map;
    //     for(int i=0;i<nx_map;i++){
    //         double x = g_xmin + i*dx_map;

    //         // obliczenie u(x,y) jako kombinacja modów
    //         double u_val = 0.0;
    //         for(int m=0;m<g_N;m++){
    //             u_val += yyyk[m] * phi[m%3](0.0,0.0); // proste przybliżenie (interpolacja)
    //         }

    //         out << u_val << " ";
    //     }
    //     out << "\n";
    // }
}

// ------------------
// PROBLEM WŁASNY: funkcja pomocnicza
// ------------------
// void solve_eigenproblem(std::vector<double>& eigvals, std::vector<std::vector<double>>& modes){
//     // prosty przykład: macierz 2x2, zastępuje rzeczywiste wyznaczanie
//     eigvals = {1.0, 2.0};
//     modes.resize(2);
//     modes[0] = std::vector<double>(g_N,1.0);
//     modes[1] = std::vector<double>(g_N,0.5);
// }