/* @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date March 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#include "../../src/cmn/cmn_error.hpp"
#include "../../src/cmn/cmn_units.hpp"
#include "../../src/phc/phc.hpp"

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

#include <boost/format.hpp>
using boost::format;

#include <blitz/array.h>
using blitz::Array;
using blitz::Range;
using blitz::firstDim;
using blitz::secondDim;

#define GNUPLOT_ENABLE_BLITZ
#include <gnuplot-iostream/gnuplot-iostream.h>

typedef double real_t;

#include "bins.hpp"

std::string zeropad(int n)
{
  std::ostringstream tmp;
  tmp << std::setw(5) << std::setfill('0') << n;
  return tmp.str();
}

int main(int argc, char **argv)
{

  string dir= argc > 1 ? argv[1] : "tmp";
 
  notice_macro("opening super-droplet netCDF file")
  NcFile nf(dir+"/out.nc", NcFile::read);

  notice_macro("opening super-droplet netCDF file")
  NcFile nf2("/users/arabas/devel/test2.nc", NcFile::read);

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

  Array<real_t, 2> rhod(nx, ny);
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

  for (size_t t = 0 ; t < nt; ++t) for (string &ext : list<string>({"eps","png"}))
  {
    notice_macro("generating frame at t=" << t)
    gp << "reset" << endl;
    // progressive-rock connoisseur palette ;)
    gp << "set palette defined (0 '#000000', 1 '#993399', 2 '#00CCFF', 3 '#66CC00', 4 '#FFFF00', 5 '#FC8727', 6 '#FD0000') maxcolors 20" << endl; 
    gp << "set view map" << endl;
    gp << "set tics out scale .25" << endl;

    gp << "set ylabel 'Y [km]'" << endl;
    gp << "set yrange [" << 0 << ":" << ny * dy/1000 << "]" << endl;
    gp << "set ytics .5" << endl;

    gp << "set contour base" << endl;
    gp << "set nosurface" << endl;
    gp << "set cntrparam levels 0" << endl;
    gp << "set nokey" << endl;

    gp << "set label 't = " << int(real_t(t) * dt_out / si::seconds) << " s' at screen .51,.99 center" << endl;
    gp << "set label 'icicle/sdm' at screen .25,.99 center" << endl;
    gp << "set label 'Zach''s 2D-bin' at screen .75,.99 center" << endl;

    int cols = 2, rows = 5;

    if (ext == "png")
      gp << "set term png enhanced size " << cols * 500 << "," << rows * 500 << endl;
    else if (ext == "eps")
      gp << "set term postscript size " << 12 * cols<< "cm," << 12 * rows << "cm solid enhanced color" << endl;
    else assert(false);

    gp << "set output '" << dir << "/test_" << zeropad(t) << "." << ext << "'" << endl;
    gp << "set multiplot layout " << rows << "," << cols << endl;

    //// spectrum plot
    {
      vector<quantity<si::length, real_t>> left_edges = bins();
      int ns = left_edges.size() - 1;

      gp << "set title 'mean particle concentration [1/cm^3]'" << endl;
      gp << "set logscale cb" << endl;
      gp << "set cbrange [1:200]" << endl;
      gp << "set xrange [0:" << ns << "]" << endl;
      gp << "set xtics (";
      for (int i = 0; i < left_edges.size()-1; i+=9)
        gp << (i!=0?",":"") << "\"" << format("%4.3g") % real_t(left_edges[i] / si::metres / real_t(1e-6)) << "\" " << i;
      gp << ")" << endl;
      gp << "set xlabel 'particle (wet) radius [{/Symbol m}m]'" << endl;

      Array<real_t, 2> tmps(ns, ny), tmp(nx, ny);

      {
        // super-droplets
        for (int i = 0; i < ns; ++i)
        {
          ostringstream name;
          name << "m_0_" << i;
          nf.getVar(name.str()).getVar(start({t,0,0,0}), count({1,nx,ny,1}), tmp.data());
          for (int y = 0; y < ny; ++y) 
            tmps(i, y) = sum(tmp(blitz::Range::all(), y)) / nx;
        }
        tmps /= 1e6; // 1/m3 -> 1/cm3
        gp << "splot '-' binary" << gp.binfmt(tmps)
          << " dx=1 dy=" << dy/1000 << " origin=(.5," << dy/2000 << ",0)"
          << " using 1 with image notitle" << endl;
        gp.sendBinary(tmps);
      }
      // 2D-bin
      { 
        for (int i = 0; i < ns; ++i)
        {
          nf2.getVar("QCDIST").getVar(start({t,0,0,i}), count({1,ny,nx,1}), tmp.data());
          for (int y = 0; y < ny; ++y) 
            tmps(i, y) = sum(tmp(y, blitz::Range::all())) / nx;
        }
        nf2.getVar("GAC").getVar(start({t,0,0}), count({1,ny,nx}), tmp.data());
        tmp.transposeSelf(secondDim, firstDim);
        tmps *= tmp; // kg-1 -> m-3
        tmps /= 1e6; // 1/m3 -> 1/cm3
        gp << "splot '-' binary" << gp.binfmt(tmps)
          << " dx=1 dy=" << dy/1000 << " origin=(.5," << dy/2000 << ",0)"
          << " using 1 with image notitle" << endl;
        gp.sendBinary(tmps);
      }
    }

    gp << "set xlabel 'X [km]'" << endl;
    gp << "set xrange [" << 0 << ":" << nx * dx/1000 << "]" << endl;
    gp << "set xtics .5" << endl;
    gp << "unset logscale cb" << endl;

    //// water vapour mixing ratio
    {
      Array<real_t, 2> tmp(nx, ny);

      gp << "set title 'water vapour mixing ratio [g/kg]'" << endl;
      gp << "set cbrange [6:8]" << endl;
      // super-droplets
      nf.getVar("rhod_rv").getVar(start({t,0,0,0}), count({1,nx,ny,1}), tmp.data()); 
      tmp /= rhod;
      gp << "splot '-' binary" << gp.binfmt(tmp) << dxdy << " using ($1*1000) with image notitle" << endl;
      gp.sendBinary(tmp);
      // 2D-bin
      nf2.getVar("QV").getVar(start({t,0,0}), count({1,ny,nx}), tmp.data()); 
      tmp.transposeSelf(secondDim, firstDim);
      gp << "splot '-' binary" << gp.binfmt(tmp) << dxdy << " using ($1*1000) with image notitle" << endl;
      gp.sendBinary(tmp);
    }


    //// potential temperature
    {
      Array<real_t, 2> tmp(nx, ny);

      gp << "set title 'potential temperature [K]'" << endl;
      gp << "set cbrange [288:293]" << endl;
      // super-droplet
      nf.getVar("rhod_th").getVar(start({t,0,0,0}), count({1,nx,ny,1}), tmp.data()); 
      tmp /= rhod;
      gp << "splot '-' binary" << gp.binfmt(tmp) << dxdy << " with image notitle" << endl;
      gp.sendBinary(tmp);
      // 2D-bin
      nf2.getVar("THETA").getVar(start({t,0,0}), count({1,ny,nx}), tmp.data()); 
      tmp.transposeSelf(secondDim, firstDim);
      gp << "splot '-' binary" << gp.binfmt(tmp) << dxdy << " with image notitle" << endl;
      gp.sendBinary(tmp);
    }

 
    //// super-droplet concentration
    {
      Array<real_t, 2> tmp(nx, ny);

      gp << "set title 'super-droplet conc. [1/dx/dy/dz]'" << endl;
      gp << "set cbrange [0:512]" << endl;
      nf.getVar("sd_conc").getVar(start({t,0,0,0}), count({1,nx,ny,1}), tmp.data()); 
      gp << "splot '-' binary" << gp.binfmt(tmp) << dxdy << " using 1 with image notitle";
      gp << endl;
      gp.sendBinary(tmp);
      // 2D-bin
      gp << "splot 0" << endl;
    }


    // effective radius
    {
      Array<real_t, 2> tmp0(nx, ny), tmp1(nx, ny);

      gp << "set title 'effective radius [{/Symbol m}m] (FSSP range)'" << endl;
      gp << "set cbrange [1:30]" << endl;
      nf.getVar("m_3").getVar(start({t,0,0,0}), count({1,nx,ny,1}), tmp0.data()); 
      nf.getVar("m_2").getVar(start({t,0,0,0}), count({1,nx,ny,1}), tmp1.data()); 
      tmp0 /= tmp1;
      tmp0 *= 1e6;
      gp << "splot '-' binary" << gp.binfmt(tmp0) << dxdy << " using 1 with image notitle";
      gp << endl;
      gp.sendBinary(tmp0);
      // 2D-bin
      gp << "splot 0" << endl;
    }

    gp << "unset multiplot" << endl;
    gp << "unset label" << endl;
  }

  string cmd="convert -monitor -delay 10 -loop 1 " + dir + "/test_*.png " + dir + "/todo.gif 1>&2";
  int status = system(cmd.c_str());
  notice_macro("done.")
}
