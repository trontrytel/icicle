#include "common.hpp"
#include "gnuplot.hpp"
#include "hdf5.hpp"

#include <map>

string zeropad(int n)
{
  ostringstream tmp;
  tmp << std::setw(3) << std::setfill('0') << n;
  return tmp.str();
}

int main(int ac, char** av)
{
  if (ac != 2) error_macro("expecting 1 argument: CMAKE_BINARY_DIR")

  std::string
    dir = string(av[1]) + "/tests/fig_a/",
    h5  = dir + "out_lgrngn.h5",
    svg = dir + "out_lgrngn_spec.svg";

  // focus to the gridbox from where the size distribution is plotted
  int focus_x = 10, focus_y = 10;
  std::map<float, float> focus_d;
  std::map<float, float> focus_w;

  //info on the number and location of histogram edges
  vector<quantity<si::length>> left_edges_rd = bins_dry();
  int nsd = left_edges_rd.size() - 1;
  vector<quantity<si::length>> left_edges_rw = bins_wet();
  int nsw = left_edges_rw.size() - 1;

  for (int i = 0; i < nsd; ++i)
  {
    const string name = "rd_rng" + zeropad(i) + "_mom0";
    blitz::Array<float, 2> tmp_d(1e-6 * h5load(h5, name));

    std::cerr << "focus_d[" << left_edges_rd[i] / 1e-6 / si::metres << "] = " << sum(tmp_d(
      blitz::Range(focus_x-1, focus_x+1),
      blitz::Range(focus_y-1, focus_y+1)
    )) / 9 << std::endl;

    focus_d[left_edges_rd[i] / 1e-6 / si::metres] = sum(tmp_d(
      blitz::Range(focus_x-1, focus_x+1),
      blitz::Range(focus_y-1, focus_y+1)
    )) / 9;
  }

  for (int i = 0; i < nsw; ++i)
  {
    const string name = "rw_rng" + zeropad(i) + "_mom0";
    blitz::Array<float, 2> tmp_w(1e-6 * h5load(h5, name));

    focus_w[left_edges_rw[i] / 1e-6 / si::metres] = sum(tmp_w(
      blitz::Range(focus_x-1, focus_x+1),
      blitz::Range(focus_y-1, focus_y+1)
    )) / 9;
  }

  notice_macro("setting-up plot parameters");
  Gnuplot gp;

  gp << "set term svg dynamic enhanced fsize 13 size 500 , 500 \n";
  gp << "set size square\n";
  gp << "set view map\n";
  gp << "set output '" << svg << "'\n";
  gp << "set logscale xy" << endl;
  gp << "set xrange [.001:100]" << endl;
  gp << "set xlabel 'particle radius [um]'" << endl;
  gp << "set yrange [.001:200]" << endl;
  gp << "set ylabel 'particle concentration (per bin) [1/cm^3]'" << endl;
  gp << "set grid" << endl;
  gp << "plot"
     << "'-' with steps title 'wet radius' lw 2 lc rgb 'blue',"
     << "'-' with steps title 'dry radius' lw 2 lc rgb 'red'" << endl;
  gp.send(focus_w);
  gp.send(focus_d);
}
