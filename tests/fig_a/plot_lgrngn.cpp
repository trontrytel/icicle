#include "common.hpp"
#include "gnuplot.hpp"
#include "hdf5.hpp"

int main(int ac, char** av)
{
  if (ac != 2) error_macro("expecting 1 argument: CMAKE_BINARY_DIR")

  std::string
    dir = string(av[1]) + "/tests/fig_a/",
    h5  = dir + "out_lgrngn.h5",
    svg = dir + "out_lgrngn.svg";

  Gnuplot gp;
  init(gp, svg, 3, 2);

  char lbl = 'a';
  for (auto &fcs : std::set<std::set<std::pair<int, int>>>({focus.first, focus.second}))
  {
    for (auto &pr : fcs) 
    {
      auto &x = pr.first;
      auto &y = pr.second;

//TODO: make the square and label more visible somehow??

      // square
      gp << "set arrow from " << x-1 << "," << y-1 << " to " << x+2 << "," << y-1 << " nohead lw 1 front\n";
      gp << "set arrow from " << x-1 << "," << y+2 << " to " << x+2 << "," << y+2 << " nohead lw 1 front\n";
      gp << "set arrow from " << x-1 << "," << y-1 << " to " << x-1 << "," << y+2 << " nohead lw 1 front\n";
      gp << "set arrow from " << x+2 << "," << y-1 << " to " << x+2 << "," << y+2 << " nohead lw 1 front\n";
      // cross 
      //gp << "set arrow from " << x-2 << "," << y+.5 << " to " << x+3 << "," << y+.5 << " nohead lw 1 front\n";
      //gp << "set arrow from " << x+.5 << "," << y-2 << " to " << x+.5 << "," << y+3 << " nohead lw 1 front\n";
      // labels
      gp << "set label " << int(lbl) << " '" << lbl << "' at " << x+(lbl%2?-8:+4) << "," << y+.5 << " front\n";

      ++lbl;
    }
  }

  // liquid water content
  { //                                                   rho_w  kg2g
    auto tmp = h5load(h5, "rw_rng000_mom3") * 4./3 * 3.14 * 1e3 * 1e3;
    gp << "set title 'liquid water content [g/m^3]'\n";
    gp << "set cbrange [0:1.5]\n";
    plot(gp, tmp);
  }

  // rain water content
  { //                                                   rho_w  kg2g
    auto tmp = h5load(h5, "rw_rng001_mom3") * 4./3 * 3.14 * 1e3 * 1e3;
    gp << "set logscale cb\n";
    gp << "set title 'rain water content [g/m^3]'\n";
    gp << "set cbrange [1e-2:1]\n";
    plot(gp, tmp);
    gp << "unset logscale cb\n";
  }

  // cloud particle concentration
  {
    auto tmp = 1e-6 * h5load(h5, "rw_rng000_mom0");
    gp << "set title 'cloud droplet concentration [cm^{-3}]'\n";
    gp << "set cbrange [0:150]\n";
    plot(gp, tmp);
  }

  // rain particle concentration
  {
    auto tmp = 1e-6 * h5load(h5, "rw_rng001_mom0");
    gp << "set title 'rain drop concentration [cm^{-3}]'\n";
    gp << "set cbrange [.01:10]\n";
    gp << "set logscale cb\n";
    plot(gp, tmp);
    gp << "unset logscale cb\n";
  }

  // effective radius
  {
    auto r_eff = h5load(h5, "rw_rng000_mom3") / h5load(h5, "rw_rng000_mom2") * 1e6;
    gp << "set title 'effective radius [um]'\n"; // TODO: Symbol nie dziala...
    gp << "set cbrange [1:20]\n";
    plot(gp, r_eff);
  }

  // radar reflectivity
  {
    auto m0 = h5load(h5, "rw_rng001_mom0");
    auto m6 = h5load(h5, "rw_rng001_mom6");
    float minval = -80, maxval = -20;
    gp << "set cbrange [" << minval << ":" << maxval << "]\n";
    auto dbZ = where(
      // if
      m0==0, 
      // then
      minval,
      // else
      10 * log10(    // reflectivity -> decibels of reflectivity
        pow(2e3,6) * // radii in metres -> diameters in milimetres
        m6 / m0      // sixth moment per unit volume
      )
    );
    gp << "set title 'radar reflectivity [dB]'\n";
    plot(gp, dbZ);
  }
}
