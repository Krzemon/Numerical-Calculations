#include "zadania.hpp"
#include "fun.hpp"
#include <vector>

void calc_zad_1(std::ofstream &out) {
    auto [nodes, elements] = generateStructuredMesh(10, 10, -5, 5, -5, 5);

    out << "# x y\n";
    for (auto &p : nodes)
        out << p.x << " " << p.y << "\n";

    out << "\n# elementy i j k\n";
    for (auto &e : elements)
        out << e.i << " " << e.j << " " << e.k << "\n";
}

void calc_zad_2(std::ofstream &out) {
    out << "# zeta eta phi0 phi1 phi2\n";

    for (double z=-1; z<=1.0001; z+=0.2)
        for (double e=-1; e<=1.0001; e+=0.2) {
            double phi0, phi1, phi2;
            shapeFunctions(z, e, phi0, phi1, phi2);
            out << z << " " << e << " " << phi0 << " " << phi1 << " " << phi2 << "\n";
        }
}

void calc_zad_3(std::ofstream &out) {
    auto [nodes, elements] = generateStructuredMesh(10, 10, -5, 5, -5, 5);

    out << "# E_ij dla 4 elementow\n";

    for (int m = 0; m < 4; m++) {
        auto &el = elements[m];
        double Ke[3][3];
        computeLocalStiffness(nodes, el, Ke);

        out << "\n# element " << m << "\n";
        for (int i=0;i<3;i++){
            for (int j=0;j<3;j++) out << Ke[i][j] << " ";
            out << "\n";
        }
    }
}

void calc_zad_4(std::ofstream &out) {
    auto [nodes, elements] = generateStructuredMesh(10, 10, -5, 5, -5, 5);

    out << "# Fj dla 4 elementow\n";

    for (int m = 0; m < 4; m++) {
        auto &el = elements[m];
        double Fe[3];
        computeLocalLoad(nodes, el, Fe);

        out << "\n# element " << m << "\n";
        out << Fe[0] << " " << Fe[1] << " " << Fe[2] << "\n";
    }
}

#include <Eigen/Dense>

void calc_zad_5(std::ofstream &out) {
    int Nx = 10, Ny = 10;
    auto [nodes, elements] = generateStructuredMesh(Nx, Ny, -5, 5, -5, 5);

    int N = nodes.size();

    Eigen::MatrixXd S = Eigen::MatrixXd::Zero(N, N);
    Eigen::VectorXd F = Eigen::VectorXd::Zero(N);

    // składanie globalnych macierzy
    for (auto &el : elements) {
        double Ke[3][3];
        double Fe[3];
        computeLocalStiffness(nodes, el, Ke);
        computeLocalLoad(nodes, el, Fe);

        int idx[3] = {el.i, el.j, el.k};
        for (int a=0;a<3;a++) {
            for (int b=0;b<3;b++) 
                S(idx[a], idx[b]) += Ke[a][b];
            F(idx[a]) += Fe[a];
        }
    }

    // warunki brzegowe Dirichleta u=0 na brzegu
    for (int i=0;i<N;i++) {
        if (isBoundaryNode(nodes[i])) {
            S.row(i).setZero();
            S.col(i).setZero();
            S(i,i)=1.0;
            F(i)=0.0;
        }
    }

    // rozwiązanie
    Eigen::VectorXd c = S.fullPivLu().solve(F);

    // zapis (x y u)
    out << "# x y u\n";
    for (int i=0;i<N;i++)
        out << nodes[i].x << " " << nodes[i].y << " " << c[i] << "\n";
}