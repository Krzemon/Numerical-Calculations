#ifndef GRID_HPP
#define GRID_HPP

#include "Types.hpp"

const double g = 9.81;

enum PointType
{
    N = 1,
    S,
    E,
    W,
    NE,
    NW,
    SE,
    SW
};

class Grid
{
public:
    Grid(double X, double Y, double dx, double dy, double dam_width, double dam_pos);
    void setInitialConditions(double h1, double h2, double u0, double v0);
    void makeTimeStep(double dt);

    void printGrid(std::ofstream &file);

private:
    struct Geometry
    {
        double _X;
        double _Y;
        double _dx;
        double _dy;
        double _dam_width;
        double _dam_pos;

        int _Nx;
        int _Ny;

        Geometry(double X, double Y, double dx, double dy, double dam_width, double dam_pos) : _X{X}, _Y{Y}, _dx(dx), _dy(dy), _dam_width{dam_width}, _dam_pos{dam_pos}
        {
            _Nx = static_cast<int>(std::round(_X / _dx)) + 1;
            _Ny = static_cast<int>(std::round(_Y / _dy)) + 1;
        }
    } _geometry;

    struct Point
    {
        VarVector _U{{0.0, 0.0, 0.0}};
        VarVector _F{{0.0, 0.0, 0.0}};
        VarVector _G{{0.0, 0.0, 0.0}};
    };

    std::vector<std::vector<Point>> _grid;
    std::vector<std::vector<int>> _state;
    std::vector<std::vector<Point>> _grid_nh_i;
    std::vector<std::vector<Point>> _grid_nh_j;

    void caclFG();
    void LaxWendroff(double dt);
};

#endif // GRID_HPP