/* @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date March 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#include <string>
using std::string;

#include <vector>
using std::vector;

#include <netcdf>
using netCDF::NcFile;
typedef vector<size_t> start;
typedef vector<size_t> count;

#include <iostream>
using std::endl;
using std::cerr;
using std::ostringstream;

#include <boost/units/systems/si.hpp>
#include <boost/units/io.hpp>
namespace si = boost::units::si;
using boost::units::quantity;

#include <blitz/array.h>
using blitz::Array;
using blitz::Range;

#define GNUPLOT_ENABLE_BLITZ
#include <gnuplot-iostream/gnuplot-iostream.h>

#define notice_macro(msg) { cerr << msg << endl; }
typedef float real_t;

std::string zeropad(int n)
{
  std::ostringstream tmp;
  tmp << std::setw(5) << std::setfill('0') << n;
  return tmp.str();
}

int main()
{
  notice_macro("opening netCDF file")
  NcFile nf("out.nc", NcFile::read);

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

  string dxdy;
  {
    ostringstream tmp;
    tmp << " dx=" << dx/1000;
    tmp << " dy=" << dy/1000;
    tmp << " origin=(" << dx/2000 << "," << dy/2000 << ",0)";
    dxdy = tmp.str();
  }

  notice_macro("allocating temp storage space")
  Array<real_t,2> tmp(nx, ny), rhod(nx, ny);
  {
    NcFile nfini("ini.nc", NcFile::read);
    for (int i = 0; i < nx; ++i)
    {
      Range j(0, ny-1);
      assert(rhod(i,j).isStorageContiguous());
      nfini.getVar("rhod").getVar(start({0,0,0}), count({1,ny,1}), rhod(i,j).data());
    }
  }

  notice_macro("setting-up plot parameters")
  Gnuplot gp;
  gp << "set term png enhanced size 1200,400" << endl;
  gp << "set view map" << endl;
  gp << "set xlabel 'X [km]'" << endl;
  gp << "set xrange [" << 0 << ":" << nx * dx/1000 << "]" << endl;
  gp << "set ylabel 'Y [km]'" << endl;
  gp << "set yrange [" << 0 << ":" << ny * dy/1000 << "]" << endl;

  gp << "set contour base" << endl;
  //gp << "set palette maxcolors 20" << endl;
  gp << "set palette defined (0 0 0 0, 1 0 0 1, 3 0 1 0, 4 1 0 0, 6 1 1 1)" << endl;
  gp << "set nosurface" << endl;
  gp << "set cntrparam levels 1" << endl;
  gp << "set nokey" << endl;

  system("mkdir -p tmp");

  for (size_t t = 0; t < nt; ++t)
  {
    notice_macro("generating frame at t=" << t)
    gp << "set label 't = " << int(real_t(t) * dt_out / si::seconds) << " s' at screen .5,.9" << endl;
    gp << "set output 'tmp/test_" << zeropad(t) << ".png'" << endl;
    gp << "set multiplot layout 1,3" << endl;

    gp << "set title 'water vapour mixing ratio [g/kg]'" << endl;
    gp << "set cbrange [5:10]" << endl;
    nf.getVar("rhod_rv").getVar(start({t,0,0,0}), count({1,nx,ny,1}), tmp.data()); 
    tmp /= rhod;
    gp << "splot '-' binary" << gp.binfmt(tmp) << dxdy << " using ($1*1000) with image notitle";
    gp << ",'-' binary" << gp.binfmt(tmp) << dxdy << " ps 0 notitle" << endl;
    gp.sendBinary(tmp);
    gp.sendBinary(tmp);

    gp << "set title 'liquid water mixing ratio [g/kg]'" << endl;
    gp << "set cbrange [0:1]" << endl;
    nf.getVar("rhod_rl").getVar(start({t,0,0,0}), count({1,nx,ny,1}), tmp.data()); 
    tmp /= rhod;
    gp << "splot '-' binary" << gp.binfmt(tmp) << dxdy << " using ($1*1000) with image notitle";
    gp << ",'-' binary" << gp.binfmt(tmp) << dxdy << " ps 0 notitle" << endl;
    gp.sendBinary(tmp);
    gp.sendBinary(tmp);

    gp << "set title 'rain water mixing ratio [g/kg]'" << endl;
    gp << "set cbrange [-.05:.05]" << endl;
    nf.getVar("rhod_rr").getVar(start({t,0,0,0}), count({1,nx,ny,1}), tmp.data()); 
    tmp /= rhod;
    gp << "splot '-' binary" << gp.binfmt(tmp) << dxdy << " using ($1*1000) with image notitle";
    gp << ",'-' binary" << gp.binfmt(tmp) << dxdy << " ps 0 notitle" << endl;
    gp.sendBinary(tmp);
    gp.sendBinary(tmp);

    gp << "unset multiplot" << endl;
    gp << "unset label 1" << endl;
  }

  system("convert -delay 10 tmp/test_*.png todo.gif");
  system("rm -rf tmp/test_*.png");
  notice_macro("done.")
}

