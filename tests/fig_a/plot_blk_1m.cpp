#include "../common.hpp"
#include "bins.hpp"
#include "gnuplot.hpp"
#include "hdf5.hpp"

int main(int ac, char** av)
{
  if (ac != 2) error_macro("expecting 1 argument: CMAKE_BINARY_DIR")

  std::string 
    dir = string(av[1]) + "/tests/fig_a/",
    h5  = dir + "out_blk_1m.h5",
    svg = dir + "out_blk_1m.svg";

  Gnuplot gp;
  init(gp, svg, 1, 2);

  {
    auto rc = h5load(h5, "rc") * 1e3;
    gp << "set title 'cloud water mixing ratio r_c [g/kg]'\n";  // TODO: -> multiply by density to compare with SD
    gp << "set cbrange [0:1.5]\n";
    plot(gp, rc);
  }

  {
    auto rr = h5load(h5, "rr") * 1e3;
    gp << "set logscale cb\n";
    gp << "set title 'rain water mixing ratio r_r [g/kg]'\n"; // TODO: -> multiply by density to compare with SD
    gp << "set cbrange [1e-2:1]\n";
    plot(gp, rr);
  }
}
