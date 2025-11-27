#pragma once
#include <vector>

struct Node { double x, y; };
struct Element { int i, j, k; };

std::pair<std::vector<Node>, std::vector<Element>>
generateStructuredMesh(int Nx, int Ny, double xmin, double xmax, double ymin, double ymax);

void shapeFunctions(double z, double e, double &phi0, double &phi1, double &phi2);

void computeLocalStiffness(const std::vector<Node> &nodes, const Element &el, double Ke[3][3]);
void computeLocalLoad(const std::vector<Node> &nodes, const Element &el, double Fe[3]);
bool isBoundaryNode(const Node &p);