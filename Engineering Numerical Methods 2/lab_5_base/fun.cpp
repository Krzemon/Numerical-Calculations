#include "fun.hpp"
#include <cmath>

// ------------------------------------------------
//  Siatka trójkątów
// ------------------------------------------------

std::pair<std::vector<Node>, std::vector<Element>>
generateStructuredMesh(int Nx, int Ny, double xmin, double xmax, double ymin, double ymax) {

    std::vector<Node> nodes;
    std::vector<Element> elements;

    double dx = (xmax - xmin) / (Nx - 1);
    double dy = (ymax - ymin) / (Ny - 1);

    // węzły
    for (int j=0;j<Ny;j++)
        for (int i=0;i<Nx;i++)
            nodes.push_back({xmin + i*dx, ymin + j*dy});

    // elementy (2 trójkąty na kwadrat)
    for (int j=0;j<Ny-1;j++)
        for (int i=0;i<Nx-1;i++) {
            int n0 = j*Nx + i;
            int n1 = n0 + 1;
            int n2 = n0 + Nx;
            int n3 = n2 + 1;

            // trójkąt 1
            elements.push_back({n0, n1, n3});
            // trójkąt 2
            elements.push_back({n0, n3, n2});
        }

    return {nodes, elements};
}

// -------------------------------------------
// Funkcje kształtu
// -------------------------------------------

void shapeFunctions(double z, double e, double &phi0, double &phi1, double &phi2) {
    phi0 = -0.5 * (e + z);
    phi1 = 0.5 * (1 + z);
    phi2 = 0.5 * (1 + e);
}

// -------------------------------------------
// Mapowanie elementu – macierz J i jej odwrotność
// -------------------------------------------

static void jacobian(const Node &A, const Node &B, const Node &C,
                     double z, double e,
                     double J[2][2], double invJ[2][2], double &detJ)
{
    double dphi_dz[3] = { -0.5, 0.5, 0.0 };
    double dphi_de[3] = { -0.5, 0.0, 0.5 };

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

// -------------------------------------------
// 7-punktowa kwadratura do trójkąta
// -------------------------------------------

static const double GP[7][2] = {
    { 1.0/3, 1.0/3 },
    { 0.2, 0.2 },
    { 0.6, 0.2 },
    { 0.2, 0.6 },
    { 0.2, 0.2 },
    { 0.6, 0.2 },
    { 0.2, 0.6 }
};
static const double GW[7] = {
    0.225,
    0.1323941527, 0.1323941527, 0.1323941527,
    0.1259391805, 0.1259391805, 0.1259391805
};

// -------------------------------------------
// Macierz sztywności elementu
// -------------------------------------------

void computeLocalStiffness(const std::vector<Node> &nodes, const Element &el, double Ke[3][3]) {
    for (int i=0;i<3;i++) for (int j=0;j<3;j++) Ke[i][j] = 0;

    Node A = nodes[el.i];
    Node B = nodes[el.j];
    Node C = nodes[el.k];

    double dphi_dz[3] = { -0.5, 0.5, 0.0 };
    double dphi_de[3] = { -0.5, 0.0, 0.5 };

    for (int k=0;k<7;k++) {

        double z = GP[k][0], e = GP[k][1], w = GW[k];

        double J[2][2], invJ[2][2], detJ;
        jacobian(A, B, C, z, e, J, invJ, detJ);

        double dphi_dx[3], dphi_dy[3];
        for (int i=0;i<3;i++) {
            dphi_dx[i] = invJ[0][0]*dphi_dz[i] + invJ[0][1]*dphi_de[i];
            dphi_dy[i] = invJ[1][0]*dphi_dz[i] + invJ[1][1]*dphi_de[i];
        }

        for (int i=0;i<3;i++)
            for (int j=0;j<3;j++)
                Ke[i][j] += w * detJ * (dphi_dx[i]*dphi_dx[j] + dphi_dy[i]*dphi_dy[j]);
    }
}

// -------------------------------------------
// Wektor obciążeń
// -------------------------------------------

static double rho(double x, double y) {
    return std::exp(-(x*x + y*y) / 2.0);
}

void computeLocalLoad(const std::vector<Node> &nodes, const Element &el, double Fe[3]) {
    for (int i=0;i<3;i++) Fe[i] = 0;

    Node A = nodes[el.i];
    Node B = nodes[el.j];
    Node C = nodes[el.k];

    for (int k=0;k<7;k++) {

        double z = GP[k][0], e = GP[k][1], w = GW[k];

        double phi0, phi1, phi2;
        shapeFunctions(z, e, phi0, phi1, phi2);
        double phi[3] = {phi0, phi1, phi2};

        double J[2][2], invJ[2][2], detJ;
        jacobian(A, B, C, z, e, J, invJ, detJ);

        double x = phi0*A.x + phi1*B.x + phi2*C.x;
        double y = phi0*A.y + phi1*B.y + phi2*C.y;

        double rh = rho(x,y);

        for (int i=0;i<3;i++)
            Fe[i] += w * detJ * phi[i] * rh;
    }
}

bool isBoundaryNode(const Node &p) {
    return (std::abs(p.x + 5) < 1e-12 ||
            std::abs(p.x - 5) < 1e-12 ||
            std::abs(p.y + 5) < 1e-12 ||
            std::abs(p.y - 5) < 1e-12);
}