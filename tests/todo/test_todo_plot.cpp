/* @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date March 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#include <vector>
#include <netcdf>
#include <iostream>

#include <boost/units/systems/si.hpp>
namespace si = boost::units::si;
using boost::units::quantity;

#include <blitz/array.h>
#define GNUPLOT_ENABLE_BLITZ
#include <gnuplot-iostream.h>

#define notice_macro(msg) { std::cerr << msg << std::endl; }
typedef float real_t;
typedef std::vector<size_t> start;
typedef std::vector<size_t> count;

template <typename T, int d>
Gnuplot &operator<<(Gnuplot &gp, const blitz::Array<T,d> &arr)
{
  gp.write(reinterpret_cast<const char*>(arr.data()), arr.size() * sizeof(T));
  return gp;
}

std::string binfmt(const blitz::Array<float,2> &arr, float dx, float dy)
{
  std::ostringstream tmp;
  tmp << " format='%float'";
  tmp << " array=(" << arr.extent(0) << "," << arr.extent(1) << ")";
  if (arr.isMajorRank(0)) tmp << "scan=yx"; // i.e. C-style ordering
  tmp << " dx=" << dx << " dy=" << dy;
  tmp << " origin=(" << dx/2 << "," << dy/2 << ",0)";
  return tmp.str();
}

int main()
{
  notice_macro("opening netCDF file")
  netCDF::NcFile nf("out.nc", netCDF::NcFile::read);

  notice_macro("reading dt_out")
  quantity<si::time, real_t> dt_out;
  {
    real_t tmp;
    nf.getVar("dt_out").getVar(&tmp);
    dt_out = tmp * si::seconds;
  }

  notice_macro("reading nx, ny")
  size_t 
    nt = nf.getDim("time").getSize(),
    nx = nf.getDim("X").getSize(),
    ny = nf.getDim("Y").getSize();

  notice_macro("reading dx, dy")
  real_t dx, dy;
  {
    real_t tmp;
    nf.getVar("X").getVar(start({1}), count({1}), &dx);
    nf.getVar("X").getVar(start({0}), count({1}), &tmp);
    dx -= tmp;
    nf.getVar("Y").getVar(start({1}), count({1}), &dy);
    nf.getVar("Y").getVar(start({0}), count({1}), &tmp);
    dy -= tmp;
  }

  notice_macro("allocating temp storage space")
  blitz::Array<real_t,2> rhod_rl(nx, ny);

  notice_macro("setting-up plot parameters")
  Gnuplot gp;
  gp << "set term png" << std::endl;
  gp << "set pm3d map" << std::endl;
  gp << "set palette" << std::endl;
  gp << "set xlabel 'X [m]'" << std::endl;
  gp << "set xrange [" << 0 << ":" << nx * dx << "]" << std::endl;
  gp << "set ylabel 'Y [m]'" << std::endl;
  gp << "set yrange [" << 0 << ":" << ny * dy << "]" << std::endl;

  system("mkdir -p tmp");

  for (size_t t = 0; t < 5; ++t)
  {
    gp << "set output 'tmp/test_" << t << ".png'" << std::endl;
    gp << "splot '-' binary" << binfmt(rhod_rl, dx, dy) << " with image notitle" << std::endl;
    nf.getVar("rhod_rl").getVar(start({t,0,0,0}), count({1,nx,ny,1}), rhod_rl.data());
    gp << rhod_rl;
  }

  system("convert tmp/test_*.png todo.gif");
  system("rm -rf tmp");

  notice_macro("done.")
}
