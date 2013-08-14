#include "common.hpp"
#include "gnuplot.hpp"
#include "hdf5.hpp"

int main(int ac, char** av)
{
  if (ac != 2) error_macro("expecting 1 argument: CMAKE_BINARY_DIR")

  string file = string(av[1]) + "/tests/fig_a/out_blk_1m.h5";

  Gnuplot gp;
  gp << "set term svg size 700, 350 fname \"Verdana\" fsize 9 \n";
  // progressive-rock connoisseur palette ;)
  gp << "set palette defined (0 '#FFFFFF', 1 '#993399', 2 '#00CCFF', 3 '#66CC00', 4 '#FFFF00', 5 '#FC8727', 6 '#FD0000')\n";
  gp << "set view map\n";

  gp << "set output '" << av[1] << "/tests/fig_a/blk_1m.svg'\n";
  gp << "set multiplot layout 1,2\n";

  gp << "set xlabel \"x/dx\"\n";
  gp << "set ylabel \"y/dy\"\n";

  auto rhod_rc = h5load(file, "rhod_rc") * 1000;
  gp << "set title \"rhod_rc [g/m3]\" \n"; 
  plot(gp, rhod_rc);
  auto rhod_rr = h5load(file, "rhod_rr") * 1000;
  gp << "set title \"rhod_rr [g/m3]\" \n"; 
  plot(gp, rhod_rr);

  gp << "set term x11\n";
}
