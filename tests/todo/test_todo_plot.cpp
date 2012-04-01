/* @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date March 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#include <vector>
#include <netcdf>

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
  //gp << "set term postscript color" << std::endl;
  gp << "set term png" << std::endl;
  gp << "set pm3d map" << std::endl;
  gp << "set palette" << std::endl;
  gp << "set xlabel 'X [m]'" << std::endl;
  gp << "set xrange [" << 0 << ":" << nx * dx << "]" << std::endl;
  gp << "set ylabel 'Y [m]'" << std::endl;
  gp << "set yrange [" << 0 << ":" << ny * dy << "]" << std::endl;

  system("mkdir tmp");

  for (size_t t = 0; t < 5; ++t)
  {
    //gp << "set output 'test_" << t << ".ps'" << std::endl;
    gp << "set output 'tmp/test_" << t << ".png'" << std::endl;
    gp << "splot '-' binary";
    gp << " format='%float' array=(" << nx << "," << ny << ") scan=yx";
    gp << " dx=" << dx << " dy=" << dy << " origin=(" << dx/2 << "," << dy/2 << ",0)";
    gp << " with image notitle" << std::endl;
    nf.getVar("rhod_rl").getVar(start({t,0,0,0}), count({1,nx,ny,1}), rhod_rl.data());
    gp.write(reinterpret_cast<const char*>(rhod_rl.data()), nx * ny * sizeof(real_t));
  }

  system("convert tmp/test_*.png todo.gif");
  system("rm -rf tmp");

  notice_macro("done.")
}
