#include "common.hpp"
#include "gnuplot.hpp"
#include "hdf5.hpp"

int main(int ac, char** av)
{
  if (ac != 2) error_macro("expecting 1 argument: CMAKE_BINARY_DIR")

  string file = string(av[1]) + "/tests/fig_a/out_lgrngn.h5";

  Gnuplot gp;

  gp << "set term svg size 700, 700 fname \"Verdana\" fsize 9 \n";
  // progressive-rock connoisseur palette ;)
  gp << "set palette defined (0 '#FFFFFF', 1 '#993399', 2 '#00CCFF', 3 '#66CC00', 4 '#FFFF00', 5 '#FC8727', 6 '#FD0000')\n";
  gp << "set view map\n";

  gp << "set output '" << av[1] << "/tests/fig_a/lgrngn.svg'\n";

  gp << "set multiplot layout 2,2\n";
  gp << "set xlabel \"x/dx\"\n";
  gp << "set ylabel \"y/dy\"\n";

  auto r_eff = h5load(file, "rw_rng0_mom3") / h5load(file, "rw_rng0_mom2") * 1e6;
  gp << "set title \"r_eff [um]\" \n";
  plot(gp, r_eff);
  auto N = h5load(file, "rw_rng0_mom0");
  gp << "set title \"N [1/m3]\" \n";
  plot(gp, N);
  auto sd_conc = h5load(file, "sd_conc");
  gp << "set title \"sd_conc [1/m3]\" \n";
  plot(gp, sd_conc);
  auto rhod_rl = h5load(file, "rw_rng0_mom3") * 4./3 * 3.14 * 1 * 1000;
  gp << "set title \"rhod_rl [g/m3]\" \n";
  plot(gp, rhod_rl);
}
