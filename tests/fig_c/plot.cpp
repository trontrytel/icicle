#include "../common.hpp"
#include "../fig_a/gnuplot.hpp"
#include "../../src/icmw8_case1.hpp"

int main(int ac, char** av)
{
  Gnuplot gp;
  init(gp, "plot.svg", 1, 1);

  using namespace icmw8_case1;

  int nx = 15, nz = 15;
  float dx = (X / si::metres) / nx, dz = (Z / si::metres) / nz;

  blitz::Array<float, 2> tmp(nx, nz);

  float A = (w_max / si::metres_per_second) * (X / si::metres) / pi<float>();

  blitz::firstIndex ix;
  blitz::secondIndex jx;

  std::vector<std::vector<float>> v(4);

  tmp = (ix + .5) / nx * (X / si::metres) / 20; 

  for (const auto &i : tmp) v[0].push_back(i);
  
  tmp = (jx + .5) / nz * (Z / si::metres) / 20;

  for (const auto &i : tmp) v[1].push_back(i);

  tmp = - A * ( 
    psi((ix+.5)/nx, (jx+.5+.5)/nz)-
    psi((ix+.5)/nx, (jx+.5-.5)/nz)
  ) / dz                // numerical derivative
  / rhod()((jx+.5) * dz); // psi defines rho_d times velocity

  for (const auto &i : tmp) v[2].push_back(i);

  tmp = A * ( 
    psi((ix+.5+.5)/nx, (jx+.5)/nz) -
    psi((ix+.5-.5)/nx, (jx+.5)/nz)
  ) / dx  
  / rhod()(jx * dz);

  for (const auto &i : tmp) v[3].push_back(i);

  gp << "scl = 3.33\n";
  gp << "set xrange [0:75]\n";
  gp << "set yrange [0:75]\n";
  gp << "set cbrange [0:5]\n";
  gp << "splot '-' using 1:2:(0):(scl*$3):(scl*$4):(0):(sqrt($3**2 + $4**2)) with vectors linecolor rgbcolor variable notitle\n";
  gp.send(v);
}
