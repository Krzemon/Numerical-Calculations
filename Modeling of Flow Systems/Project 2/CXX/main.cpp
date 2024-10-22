#include "Grid.hpp"

// strumien=wysokosc*predkosc
// dx=1m
// obszar 100 na 100
// najpierw po lewej full szerokosc 10m wysokosci
// w polowie tama 40 m szerokosci
// po prawej rzeka na 40 m szerokosci i 5m wysokosci

const int N_frames = 100;
const double dt = 0.05;

int main()
{
  double L = 100.0;
  double step = 1.0;
  double dam_width = 0.4 * L;
  double dam_pos = 0.5 * L;

  Grid grid1(L, L, step, step, dam_width, dam_pos);

  double h1 = 10.0;
  double h2 = 5.0;
  double u0 = 0.0;
  double v0 = 0.0;

  grid1.setInitialConditions(h1, h2, u0, v0);

  std::ofstream grid_file("data/grid.dat", std::ios::out);

  for (int i = 0; i < N_frames; i++)
  {
    grid1.printGrid(grid_file);
    grid1.makeTimeStep(dt);
  }
}
