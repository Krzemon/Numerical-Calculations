#include <functional>
#include <fstream>
#include <string>
#include <cmath>
#include <array>

const double PI = 3.141592653589793238463;

int main(int argc, char *argv[])
{
  const double mu0 = 1;
  const double epsilon0 = 1;
  const double R = 1;
  const double sigma = 1;
  const double omega = 0.1;
  const double alpha = PI / 4;

  const size_t N = 201;
  const size_t M = 201;

  const double dtheta = PI / N;
  const double dphi = 2 * PI / M;

  const std::array<double, 3> omega_vec = {omega * sin(alpha), omega * cos(alpha), 0};

  const double L = 3.0;
  const size_t K = 41;
  const double dx = 2 * L / K;
  const double dy = 2 * L / K;

  auto coeff = [](size_t i, size_t N)
  {
    if (i == 0 || i == N)
      return 1;
    else if (i % 2 == 0)
      return 2;
    else
      return 4;
  };

  auto cross_product = [](const std::array<double, 3> &a, const std::array<double, 3> &b)
  {
    std::array<double, 3> c;
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
    return c;
  };

  double x, y, z, V, ex, ey, bx, by;
  std::ofstream output_file("../data/output.dat", std::ios::out);

  for (size_t ix = 0; ix <= K; ix++)
  {
    for (size_t jy = 0; jy <= K; jy++)
    {

      x = -L + ix * dx;
      y = -L + jy * dy;
      z = 0;

      V = 0;

      ex = ey = 0;
      bx = by = 0;

      for (size_t i = 0; i < N; i++)
      {
        for (size_t j = 0; j < M; j++)
        {
          double theta1 = i * dtheta;
          double phi1 = j * dphi;
          double r = sqrt(x * x + y * y + z * z);
          double x1 = R * sin(theta1) * cos(phi1);
          double y1 = R * sin(theta1) * sin(phi1);
          double z1 = R * cos(theta1);
          std::array<double, 3> r1_vec = {x1, y1, z1};

          double r_norm = sqrt((x - x1) * (x - x1) + (y - y1) * (y - y1) + (z - z1) * (z - z1));
          double W = sigma * R * R / 4 / PI * dtheta / 3 * dphi / 3;
          V += W / epsilon0 * coeff(i, N) * coeff(j, M) * sin(theta1) / r_norm;
          ex += W / epsilon0 * coeff(i, N) * coeff(j, M) * sin(theta1) * (x - x1) / pow(r_norm, 3);
          ey += W / epsilon0 * coeff(i, N) * coeff(j, M) * sin(theta1) * (y - y1) / pow(r_norm, 3);
          std::array<double, 3> g1_vec = cross_product(omega_vec, r1_vec);
          std::array<double, 3> g2_vec = cross_product({x - x1, y - y1, z - z1}, g1_vec);
          bx += -mu0 * W * coeff(i, N) * coeff(j, M) * sin(theta1) * g2_vec[0] / pow(r_norm, 3);
          by += -mu0 * W * coeff(i, N) * coeff(j, M) * sin(theta1) * g2_vec[1] / pow(r_norm, 3);
        }
      }

      output_file << x << " " << y << " " << V << " " << ex << " " << ey << " " << bx << " " << by << std::endl;
    }
    
    output_file << std::endl;
  }

  output_file.close();
}