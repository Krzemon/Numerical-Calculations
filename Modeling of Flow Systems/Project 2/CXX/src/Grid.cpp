#include "Grid.hpp"

#include <iostream>

Grid::Grid(double X, double Y, double dx, double dy, double dam_width, double dam_pos) : _geometry(Geometry(X, Y, dx, dy, dam_width, dam_pos))
{
    _grid.resize(_geometry._Nx, std::vector<Point>(_geometry._Ny, Point()));
    _state.resize(_geometry._Nx, std::vector<int>(_geometry._Ny, -1));
    _grid_nh_i.resize(_geometry._Nx - 1, std::vector<Point>(_geometry._Ny, Point()));
    _grid_nh_j.resize(_geometry._Nx, std::vector<Point>(_geometry._Ny - 1, Point()));
}

void Grid::setInitialConditions(double h1, double h2, double u0, double v0)
{
    int dam_i = static_cast<int>(_geometry._dam_pos / _geometry._dx);
    int river_bottom_j = static_cast<int>((_geometry._Y - _geometry._dam_width) / 2 / _geometry._dy);
    int river_top_j = static_cast<int>((_geometry._Y + _geometry._dam_width) / 2 / _geometry._dy);

#pragma omp parallel for collapse(2)
    for (int i = 0; i < _geometry._Nx; i++)
    {
        for (int j = 0; j < _geometry._Ny; j++)
        {

            if (i == 0)
            {
                _grid[i][j]._U = {h1, u0, v0};
                _state[i][j] = W;
            }
            else if (i == dam_i)
            {
                _grid[i][j]._U = {h1, u0, v0};
                if (j >= river_bottom_j && j <= river_top_j)
                {
                    _state[i][j] = 0;
                }
                else
                {
                    _state[i][j] = E;
                }
            }
            else if (i == _geometry._Nx - 1)
            {

                if (j >= river_bottom_j && j <= river_top_j)
                {
                    _state[i][j] = E;
                    _grid[i][j]._U = {h2, u0, v0};
                }
                else
                {
                    _state[i][j] = -1;
                    _grid[i][j]._U = {0.0, u0, v0};
                }
            }
            else if (i < dam_i)
            {
                _grid[i][j]._U = {h1, u0, v0};
                if (j == 0)
                {
                    _state[i][j] = S;
                }
                else if (j == _geometry._Ny - 1)
                {
                    _state[i][j] = N;
                }
                else
                {
                    _state[i][j] = 0;
                }
            }
            else
            {
                if (j == river_bottom_j)
                {
                    _state[i][j] = S;
                    _grid[i][j]._U = {h2, u0, v0};
                }
                else if (j == river_top_j)
                {
                    _state[i][j] = N;
                    _grid[i][j]._U = {h2, u0, v0};
                }
                else if (j > river_bottom_j && j < river_top_j)
                {
                    _state[i][j] = 0;
                    _grid[i][j]._U = {h2, u0, v0};
                }
                else
                {
                    _state[i][j] = -1;
                    _grid[i][j]._U = {0.0, u0, v0};
                }
            }
        }
    }

    _state[0][0] = SW;
    _state[0][_geometry._Ny - 1] = NW;
    _state[dam_i][0] = SE;
    _state[dam_i][river_bottom_j] = SE;
    _state[dam_i][river_top_j] = NE;
    _state[dam_i][_geometry._Ny - 1] = NE;
    _state[_geometry._Nx - 1][river_bottom_j] = SE;
    _state[_geometry._Nx - 1][river_top_j] = NE;

    caclFG();
}

void Grid::makeTimeStep(double dt)
{
    LaxWendroff(dt);

#pragma omp parallel for collapse(2)
    for (int i = 1; i < _geometry._Nx - 1; i++)
    {
        for (int j = 1; j < _geometry._Ny - 1; j++)
        {
            if (_state[i][j] == 0)
                _grid[i][j]._U = _grid[i][j]._U - dt / _geometry._dx * (_grid_nh_i[i][j]._F - _grid_nh_i[i - 1][j]._F) - dt / _geometry._dy * (_grid_nh_j[i][j]._G - _grid_nh_j[i][j - 1]._G);
        }
    }

#pragma omp parallel for collapse(2)
    for (int i = 0; i < _geometry._Nx; i++)
    {
        for (int j = 0; j < _geometry._Ny; j++)
        {
            switch (_state[i][j])
            {
            case N:
                _grid[i][j]._U = _grid[i][j - 1]._U;
                break;
            case S:
                _grid[i][j]._U = _grid[i][j + 1]._U;
                break;
            case E:
                _grid[i][j]._U = _grid[i - 1][j]._U;
                break;
            case W:
                _grid[i][j]._U = _grid[i + 1][j]._U;
                break;
            case NE:
                _grid[i][j]._U = _grid[i - 1][j - 1]._U;
                break;
            case NW:
                _grid[i][j]._U = _grid[i + 1][j - 1]._U;
                break;
            case SE:
                _grid[i][j]._U = _grid[i - 1][j + 1]._U;
                break;
            case SW:
                _grid[i][j]._U = _grid[i + 1][j + 1]._U;
                break;
            }
        }
    }

    caclFG();
}

void Grid::printGrid(std::ofstream &file)
{
    for (int i = 0; i < _geometry._Nx; i++)
    {
        for (int j = 0; j < _geometry._Ny; j++)
        {
            double h = _grid[i][j]._U[0];
            double uh = _grid[i][j]._U[1];
            double vh = _grid[i][j]._U[2];

            file << i * _geometry._dx << '\t' << j * _geometry._dy << '\t' << _state[i][j] << '\t' << h << '\t' << ((h < 1e-6) ? 0.0 : (uh / h)) << '\t' << ((h < 1e-6) ? 0.0 : (vh / h)) << std::endl;
        }

        file << std::endl;
    }
    file << std::endl;
}

// private

void Grid::caclFG()
{
#pragma omp parallel for collapse(2)
    for (int i = 0; i < _geometry._Nx; i++)
    {
        for (int j = 0; j < _geometry._Ny; j++)
        {
            if (_state[i][j] >= 0)
            {
                double h = _grid[i][j]._U[0];
                double uh = _grid[i][j]._U[1];
                double vh = _grid[i][j]._U[2];

                _grid[i][j]._F = {uh, uh * uh / h + 0.5 * g * h * h, uh * vh / h};
                _grid[i][j]._G = {vh, uh * vh / h, vh * vh / h + 0.5 * g * h * h};
            }
        }
    }
}

void Grid::LaxWendroff(double dt)
{
#pragma omp parallel for collapse(2)
    for (int i = 0; i < _geometry._Nx; i++)
    {
        for (int j = 0; j < _geometry._Nx; j++)
        {
            if (_state[i][j] >= 0)
            {
                if (i < _geometry._Nx - 1)
                {
                    _grid_nh_i[i][j]._U = 0.5 * (_grid[i + 1][j]._U + _grid[i][j]._U) - 0.5 * dt / _geometry._dx * (_grid[i + 1][j]._F - _grid[i][j]._F);

                    double h = _grid_nh_i[i][j]._U[0];
                    double uh = _grid_nh_i[i][j]._U[1];
                    double vh = _grid_nh_i[i][j]._U[2];

                    _grid_nh_i[i][j]._F = {uh, uh * uh / h + 0.5 * g * h * h, uh * vh / h};
                }

                if (j < _geometry._Ny - 1)
                {
                    _grid_nh_j[i][j]._U = 0.5 * (_grid[i][j + 1]._U + _grid[i][j]._U) - 0.5 * dt / _geometry._dy * (_grid[i][j + 1]._G - _grid[i][j]._G);

                    double h = _grid_nh_j[i][j]._U[0];
                    double uh = _grid_nh_j[i][j]._U[1];
                    double vh = _grid_nh_j[i][j]._U[2];

                    _grid_nh_j[i][j]._G = {vh, uh * vh / h, vh * vh / h + 0.5 * g * h * h};
                }
            }
        }
    }
}