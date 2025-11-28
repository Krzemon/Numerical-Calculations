#include "solution.hpp"

void calc_zad_1(std::ofstream &out) {
    make_grid();
    out << "# index x y nx=" << nx << " ny=" << ny << "\n";
    for (const auto &p : nodes)
        out << p.index << " " << p.x << " " << p.y << "\n";
    out << "\n# elementy i j k\n";
    for (const auto &e : elements)
        out << e.index << " " << e.i << " " << e.j << " " << e.k << "\n";
    std::cout << "Wygenerowano siatke wezlow i elementow trojkatnych\n";
}

void calc_zad_2(std::ofstream &out) {
    out << std::fixed << std::setprecision(1);
    out << "# zeta eta phi0 phi1 phi2\n";
    for (double z=-1; z<=1.0001; z+=0.2)
        for (double e=-1; e<=1.0001; e+=0.2) {
            double phi0, phi1, phi2;
            shape_functions(z, e, phi0, phi1, phi2);
            out << std::setw(4) << z << " " << std::setw(4) << e << " " << std::setw(4) << phi0 << " " << std::setw(4) << phi1 << " " << std::setw(4) << phi2 << "\n";
        }
    std::cout << "Wygenerowano wartosci liniowych funckji ksztaltu w przestrzeni referencyjnej\n";
}

void calc_zad_3a(std::ofstream &out) {
    std::cout << std::fixed << std::setprecision(1);
    out << std::fixed << std::setprecision(1);
    for (auto m: {0,10,50,80}) { 
        auto &el = elements[m];
        double E_ijm[3][3];
        assemble_local_matrice_triangle(nodes, el, E_ijm);
        std::cout << "# element " << m+1 << "\n";
        out << "# element " << m+1 << "\n";
        for (int i=0; i<=2; ++i){
            for (int j=0; j<=2; ++j){
                std::cout  << std::setw(4) << E_ijm[i][j] << " ";
                out << std::setw(4) << E_ijm[i][j] << " ";
            }
            std::cout << "\n";
            out << "\n";
        }
        std::cout << "\n";
        out << "\n";
    }
}

void calc_zad_3b(std::ofstream &out) {
    std::cout << std::fixed << std::setprecision(1);
    out << std::fixed << std::setprecision(1);
    for (auto m: {0,10,50,80}) { 
        auto &el = elements[m];
        double F_jm[3];
        assemble_local_vector_triangle(nodes, el, F_jm);
        std::cout << "# element " << m+1 << "\n";
        out << "# element " << m+1 << "\n";
        for (int j=0; j<=2; ++j){
            std::cout  << std::setw(4) << F_jm[j] << " ";
            out << std::setw(4) << F_jm[j] << " ";
        }
        std::cout << "\n";
        out << "\n";
    }
}

void calc_zad_4(std::ofstream &out) {
    // int N = nx*ny;

    // // --- GLOBALNA MACIERZ I WEKTOR GSL ---
    // gsl_matrix *S = gsl_matrix_calloc(N, N);
    // gsl_vector *F = gsl_vector_calloc(N);

    // // składanie globalnych macierzy
    // for (auto &el : elements) {
    //     double Ke[3][3];
    //     double Fe[3];
    //     computeLocalStiffness(nodes, el, Ke);
    //     computeLocalLoad(nodes, el, Fe);

    //     int idx[3] = {el.i, el.j, el.k};
    //     for (int a=0; a<3; a++) {
    //         for (int b=0; b<3; b++) {
    //             double old = gsl_matrix_get(S, idx[a], idx[b]);
    //             gsl_matrix_set(S, idx[a], idx[b], old + Ke[a][b]);
    //         }
    //         double oldF = gsl_vector_get(F, idx[a]);
    //         gsl_vector_set(F, idx[a], oldF + Fe[a]);
    //     }
    // }

    // // --- WARUNKI BRZEGOWE DIRICHLETA u=0 ---
    // for (int i=0; i<N; i++) {
    //     if (isBoundaryNode(nodes[i])) {
    //         // zerujemy wiersz i kolumnę
    //         for (int j=0; j<N; j++) {
    //             gsl_matrix_set(S, i, j, 0.0);
    //             gsl_matrix_set(S, j, i, 0.0);
    //         }
    //         gsl_matrix_set(S, i, i, 1.0);
    //         gsl_vector_set(F, i, 0.0);
    //     }
    // }

    // // --- ROZWIĄZANIE S·c = F przy pomocy GSL ---
    // gsl_permutation *perm = gsl_permutation_alloc(N);
    // int signum;

    // gsl_linalg_LU_decomp(S, perm, &signum);

    // gsl_vector *c = gsl_vector_alloc(N);
    // gsl_linalg_LU_solve(S, perm, F, c);

    // // --- ZAPIS DO PLIKU ---
    // out << "# x y u\n";
    // for (int i=0; i<N; i++)
    //     out << nodes[i].x << " " << nodes[i].y << " "
    //         << gsl_vector_get(c, i) << "\n";

    // // --- CZYSZCZENIE ---
    // gsl_matrix_free(S);
    // gsl_vector_free(F);
    // gsl_vector_free(c);
    // gsl_permutation_free(perm);
}

void calc_zad_5(std::ofstream &out) {}
void calc_zad_6(std::ofstream &out) {}