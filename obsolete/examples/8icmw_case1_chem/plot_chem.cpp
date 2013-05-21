/* @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date March 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#include <string>
using std::string;

#include <list>
using std::list;

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
using boost::units::divide_typeof_helper; // TODO: to powinno by� te� inkludowane w phc.hpp
using boost::units::pow;

#include <blitz/array.h>
using blitz::Array;
using blitz::Range;

#define GNUPLOT_ENABLE_BLITZ
#include <gnuplot-iostream/gnuplot-iostream.h>

// mkdir
#include <sys/stat.h>
#include <sys/types.h>

#include "../../src/cmn/cmn_error.hpp"
#include "../../src/cmn/cmn_units.hpp"
#include "../../src/phc/phc.hpp"

typedef float real_t;

std::string zeropad(int n)
{
  std::ostringstream tmp;
  tmp << std::setw(5) << std::setfill('0') << n;
  return tmp.str();
}

int main(int argc, char **argv)
{

  string dir= argc > 1 ? argv[1] : "tmp";
 
  notice_macro("opening netCDF file")
  NcFile nf(dir+"/out.nc", NcFile::read);

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
  Array<real_t,2> tmp0(nx, ny), tmp1(nx, ny), rhod(nx, ny), rv(nx,ny), th(nx,ny);
  {
    NcFile nfini(dir+"/ini.nc", NcFile::read);
    for (int i = 0; i < nx; ++i)
    {
      Range j(0, ny-1);
      assert(rhod(i,j).isStorageContiguous());
      nfini.getVar("rhod").getVar(start({0,0,0}), count({1,ny,1}), rhod(i,j).data());
    }
  }


  notice_macro("setting-up plot parameters")
  Gnuplot gp;

  for (size_t t = 0 ; t < nt; ++t) for (string &ext : list<string>({/*"eps",*/"png"}))
  {
    notice_macro("generating frame at t=" << t)
    gp << "reset" << endl;
    // progressive-rock connoisseur palette ;)
    gp << "set palette defined (0 '#ffffff', 1 '#993399', 2 '#00CCFF', 3 '#66CC00', 4 '#FFFF00', 5 '#FC8727', 6 '#FD0000')" << endl; 
    gp << "set view map" << endl;
    gp << "set xlabel 'X [km]'" << endl;
    gp << "set xrange [" << 0 << ":" << nx * dx/1000 << "]" << endl;
    gp << "set ylabel 'Y [km]'" << endl;
    gp << "set yrange [" << 0 << ":" << ny * dy/1000 << "]" << endl;

    gp << "set contour base" << endl;
    gp << "set nosurface" << endl;
    gp << "set cntrparam levels 0" << endl;
//    gp << "set nokey" << endl;
//    gp << "set label 't = " << int(real_t(t) * dt_out / si::seconds) << " s' at screen .48,.96 left" << endl;

    string micro, cmdline;
    nf.getAtt("command_line").getValues(cmdline); 
    if (cmdline.find("--eqs todo_sdm") != string::npos)
      micro = "sdm";
    else 
      assert(false);

    if (ext == "png")
      gp << "set term png enhanced size 1200," << 700 << endl;
    else if (ext == "eps")
      gp << "set term postscript size 40cm, 25cm solid enhanced color font 'Helvetica, 12'" << endl;
    else assert(false);

    gp << "set output '" << dir << "/test_" << zeropad(t) << "." << ext << "'" << endl;
    gp << "set multiplot layout 2,3" << endl;

    gp << "set title 'water vapour mixing ratio [g/kg]'" << endl;
    gp << "set cbrange [6:8]" << endl;
    nf.getVar("rhod_rv").getVar(start({t,0,0,0}), count({1,nx,ny,1}), rv.data()); 
    rv /= rhod;
    gp << "splot '-' binary" << gp.binfmt(rv) << dxdy << " using ($1*1000) with image notitle";
    gp << endl;
    gp.sendBinary(rv);

    gp << "set title 'potential temperature [K]'" << endl;
    gp << "set cbrange [288:293]" << endl;
    nf.getVar("rhod_th").getVar(start({t,0,0,0}), count({1,nx,ny,1}), th.data()); 
    th /= rhod;
    gp << "splot '-' binary" << gp.binfmt(th) << dxdy << " with image notitle";
    gp << endl;
    gp.sendBinary(th);

//    gp << "set title 'super-droplet conc. [1/dx/dy/dz]'" << endl;
//    gp << "set cbrange [0:512]" << endl;
//    nf.getVar("sd_conc").getVar(start({t,0,0,0}), count({1,nx,ny,1}), tmp0.data()); 
//    gp << "splot '-' binary" << gp.binfmt(tmp0) << dxdy << " using 1 with image notitle";
//    gp << endl;
//    gp.sendBinary(tmp0);


    gp << "set title 'effective radius [{/Symbol m}m]'" << endl;
    gp << "set cbrange [1:30]" << endl;
    nf.getVar("m_3").getVar(start({t,0,0,0}), count({1,nx,ny,1}), tmp0.data()); 
    nf.getVar("m_2").getVar(start({t,0,0,0}), count({1,nx,ny,1}), tmp1.data()); 
    tmp0 /= tmp1;
    tmp0 *= 1e6;
    gp << "splot '-' binary" << gp.binfmt(tmp0) << dxdy << " using 1 with image notitle";
    gp << endl;
    gp.sendBinary(tmp0);

    // sdm_chem relevant plots:
    gp << "set format cb \"%g\"" << endl;

/*      gp << "set title 'mean H mass in the droplet'" << endl;
//      gp << "set autoscale cb" << endl;
      gp << "set cbrange [0:1*1e-19]" << endl;
      nf.getVar("c_H").getVar(start({t,0,0,0}), count({1,nx,ny,1}), tmp0.data()); 
      nf.getVar("n_tot").getVar(start({t,0,0,0}), count({1,nx,ny,1}), tmp1.data());     
      tmp0 = tmp0 / tmp1;
      gp << "splot '-' binary" << gp.binfmt(tmp0) << dxdy << " using 1 with image notitle";
      gp << endl;
      gp.sendBinary(tmp0);
*/
      gp << "set title 'mean SO_2 mass in the droplet'" << endl;
//      gp << "set autoscale cb" << endl;
      gp << "set cbrange [0:1*1e-22]" << endl;
      nf.getVar("c_SO2").getVar(start({t,0,0,0}), count({1,nx,ny,1}), tmp0.data()); 
      nf.getVar("n_tot").getVar(start({t,0,0,0}), count({1,nx,ny,1}), tmp1.data());     
      tmp0 = tmp0 / tmp1;
      gp << "splot '-' binary" << gp.binfmt(tmp0) << dxdy << " using 1 with image notitle";
      gp << endl;
      gp.sendBinary(tmp0);

      gp << "set title 'mean HSO_3 mass in the droplet'" << endl;
//      gp << "set autoscale cb" << endl;
      gp << "set cbrange [0:1*1e-18]" << endl;
      nf.getVar("c_HSO3").getVar(start({t,0,0,0}), count({1,nx,ny,1}), tmp0.data()); 
      nf.getVar("n_tot").getVar(start({t,0,0,0}), count({1,nx,ny,1}), tmp1.data());     
      tmp0 = tmp0 / tmp1;
      gp << "splot '-' binary" << gp.binfmt(tmp0) << dxdy << " using 1 with image notitle";
      gp << endl;
      gp.sendBinary(tmp0);

      gp << "set title 'mean S(VI) mass in the droplet'" << endl;
//      gp << "set autoscale cb" << endl;
      gp << "set cbrange [0:5*1e-19]" << endl;
      nf.getVar("c_S_VI").getVar(start({t,0,0,0}), count({1,nx,ny,1}), tmp0.data()); 
      nf.getVar("n_tot").getVar(start({t,0,0,0}), count({1,nx,ny,1}), tmp1.data());     
      tmp0 = tmp0 / tmp1;
      gp << "splot '-' binary" << gp.binfmt(tmp0) << dxdy << " using 1 with image notitle";
      gp << endl;
      gp.sendBinary(tmp0);
/*
      gp << "set title 'mean SO_4 mass in the droplet'" << endl;
//      gp << "set autoscale cb" << endl;
      gp << "set cbrange [0:5*1e-19]" << endl;
      nf.getVar("c_SO4").getVar(start({t,0,0,0}), count({1,nx,ny,1}), tmp0.data()); 
      nf.getVar("n_tot").getVar(start({t,0,0,0}), count({1,nx,ny,1}), tmp1.data());     
      tmp0 = tmp0 / tmp1;
      gp << "splot '-' binary" << gp.binfmt(tmp0) << dxdy << " using 1 with image notitle";
      gp << endl;
      gp.sendBinary(tmp0);
*/
    gp << "unset multiplot" << endl;
    gp << "unset label" << endl;
  }

//  string cmd="convert -monitor -delay 10 -loop 1 " + dir + "/test_*.png " + dir + "/todo.gif 1>&2";
//  int status = system(cmd.c_str());
//  notice_macro("done.")
}
