#include "../common.hpp"
#include "bins.hpp"
#include "gnuplot.hpp"
#include "hdf5.hpp"

int main(int ac, char** av)
{
  if (ac != 2) error_macro("expecting 1 argument: CMAKE_BINARY_DIR")

  std::string
    dir = string(av[1]) + "/tests/fig_a/",
    h5  = dir + "out_blk_2m.h5",
    svg = dir + "out_blk_2m.svg";

  Gnuplot gp;
  init(gp, svg, 2, 2);

  {
    auto rc = h5load(h5, "rc") * 1e3;
    gp << "set title 'cloud water mixing ratio r_c [g/kg]'\n"; // TODO: *rho_d
    gp << "set cbrange [0:1.5]\n";
    plot(gp, rc);
  }

  {
    auto rr = h5load(h5, "rr") * 1e3;
    gp << "set logscale cb\n";
    gp << "set title 'rain water mixing ratio r_r [g/kg]'\n"; // TODO: *rho_d
    gp << "set cbrange [1e-2:1]\n";
    plot(gp, rr);
    gp << "unset logscale cb\n";
  }

  {
    auto nc = h5load(h5, "nc") * 1e-6;
    gp << "set title 'cloud droplet specific concentration n_c [mg^{-1}]'\n"; // TODO: *rho_d
    gp << "set cbrange [0:150]\n";
    plot(gp, nc);
  }

  {
    auto nr = h5load(h5, "nr") * 1e-6;
    gp << "set title 'rain drop specific concentration n_r [mg^{-1}]'\n"; // TODO: *rho_d
    gp << "set cbrange [0.01:10]\n";
    gp << "set logscale cb\n";
    plot(gp, nr);
  }
}
