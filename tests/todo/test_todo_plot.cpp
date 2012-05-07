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
#include <boost/units/cmath.hpp>
#include <boost/units/io.hpp>
namespace si = boost::units::si;
using boost::units::quantity;
using boost::units::divide_typeof_helper; // TODO: to powinno byæ te¿ inkludowane w phc.hpp
using boost::units::pow;

#include <blitz/array.h>
using blitz::Array;
using blitz::Range;

#define GNUPLOT_ENABLE_BLITZ
#include <gnuplot-iostream/gnuplot-iostream.h>

#include "../../src/phc.hpp"

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
  Array<real_t,2> tmp(nx, ny), rhod(nx, ny), rv(nx,ny), th(nx,ny);
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
  gp << "set term png enhanced size 800,800" << endl;
  gp << "set view map" << endl;
  gp << "set xlabel 'X [km]'" << endl;
  gp << "set xrange [" << 0 << ":" << nx * dx/1000 << "]" << endl;
  gp << "set ylabel 'Y [km]'" << endl;
  gp << "set yrange [" << 0 << ":" << ny * dy/1000 << "]" << endl;

  gp << "set contour base" << endl;
  gp << "set nosurface" << endl;
  gp << "set cntrparam levels 0" << endl;
  gp << "set nokey" << endl;

  // progressive-rock connoisseur palette ;)
  gp << "set palette defined (0 '#000000', 1 '#993399', 2 '#00CCFF', 3 '#66CC00', 4 '#FFFF00', 5 '#FC8727', 6 '#FD0000')" << endl; 

  system("mkdir -p tmp");

  for (size_t t = 0; t < nt; ++t) 
  {
    notice_macro("generating frame at t=" << t)
    gp << "set label 't = " << int(real_t(t) * dt_out / si::seconds) << " s' at screen .48,.96 left" << endl;
    gp << "set output 'tmp/test_" << zeropad(t) << ".png'" << endl;
    gp << "set multiplot layout 2,2" << endl;

    gp << "set title 'water vapour mixing ratio [g/kg]'" << endl;
    gp << "set cbrange [6:8]" << endl;
    nf.getVar("rhod_rv").getVar(start({t,0,0,0}), count({1,nx,ny,1}), rv.data()); 
    rv /= rhod;
    gp << "splot '-' binary" << gp.binfmt(rv) << dxdy << " using ($1*1000) with image notitle";
    //gp << ",'-' binary" << gp.binfmt(rv) << dxdy << " ps 0 notitle";
    gp << endl;
    gp.sendBinary(rv);
    //gp.sendBinary(rv);

/*
    gp << "set title 'liquid water mixing ratio [g/kg]'" << endl;
    gp << "set cbrange [0:1]" << endl;
    nf.getVar("rhod_rl").getVar(start({t,0,0,0}), count({1,nx,ny,1}), tmp.data()); 
    tmp /= rhod;
    gp << "splot '-' binary" << gp.binfmt(tmp) << dxdy << " using ($1*1000) with image notitle";
    //gp << ",'-' binary" << gp.binfmt(tmp) << dxdy << " ps 0 notitle";
    gp << endl;
    gp.sendBinary(tmp);
    //gp.sendBinary(tmp);

    gp << "set title 'rain water mixing ratio [g/kg]'" << endl;
    gp << "set cbrange [0.:.02]" << endl;
    nf.getVar("rhod_rr").getVar(start({t,0,0,0}), count({1,nx,ny,1}), tmp.data()); 
    tmp /= rhod;
    gp << "splot '-' binary" << gp.binfmt(tmp) << dxdy << " using ($1*1000) with image notitle";
    //gp << ",'-' binary" << gp.binfmt(tmp) << dxdy << " ps 0 notitle";
    gp << endl;
    gp.sendBinary(tmp);
    //gp.sendBinary(tmp);

    gp << "set title 'potential temperature [K]'" << endl;
    gp << "set cbrange [288:293]" << endl;
    nf.getVar("rhod_th").getVar(start({t,0,0,0}), count({1,nx,ny,1}), th.data()); 
    th /= rhod;
    gp << "splot '-' binary" << gp.binfmt(th) << dxdy << " with image notitle";
    //gp << ",'-' binary" << gp.binfmt(th) << dxdy << " ps 0 notitle";
    gp << endl;
    gp.sendBinary(th);
    //gp.sendBinary(th);

    gp << "set label 'results obtained with icicle - a GPL-ed C++ MPDATA-based solver from University of Warsaw' at screen .98,.02 right" << endl;
    gp << "set label '8th International Cloud Modeling Workshop 2012: Case 1' at screen .02,.02 left" << endl;
*/

    gp << "set title 'super-droplet conc. [1/dx/dy]'" << endl;
    gp << "set cbrange [0:50]" << endl;
    nf.getVar("sd_conc").getVar(start({t,0,0,0}), count({1,nx,ny,1}), tmp.data()); 
    gp << "splot '-' binary" << gp.binfmt(tmp) << dxdy << " using 1 with image notitle";
    //gp << ",'-' binary" << gp.binfmt(tmp) << dxdy << " ps 0 notitle";
    gp << endl;
    gp.sendBinary(tmp);
    //gp.sendBinary(tmp);

    gp << "unset label" << endl;
    gp << "unset multiplot" << endl;
  }

  system("convert -delay 10 tmp/test_*.png todo.gif");
  system("rm -rf tmp/test_*.png");
  notice_macro("done.")
}

