#include "common.hpp"
#include "gnuplot.hpp"
#include "hdf5.hpp"

int main(int ac, char** av)
{
  if (ac != 2) error_macro("expecting 1 argument: CMAKE_BINARY_DIR")

  string file = string(av[1]) + "/tests/fig_a/out_blk_2m.h5";

  Gnuplot gp;
  gp << "set term svg\n";
  gp << "set view map\n";
  gp << "set multiplot layout 1,3\n";

  gp << "set output '" << av[1] << "/tests/fig_a/plot.svg'\n";

  gp << "set cbrange [.006:.01]\n";
  auto th = h5load(file, "rhod_th") / h5load(file, "rhod");
  plot(gp, th, "theta [K]");
}
