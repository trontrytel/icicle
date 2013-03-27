/* @file
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date March 2013
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

// mkdir
#include <sys/stat.h>
#include <sys/types.h>

#include "../../src/cmn/cmn_error.hpp"
#include "../../src/cmn/cmn_units.hpp"
#include "../../src/phc/phc.hpp"

typedef float real_t;

int main()
{
  // focus to the gridbox from where the size distribution is plotted
  int focus_x = 15, focus_y = 55, focus_t=10;
  std::map<real_t, real_t> focus;

  //info on the number and location of histogram edges
  //TODO get those from the command line or calc_chem
  int nsd = 30;
  int nsw = 80;
  float left_edges_rd[nsd], left_edges_rw[nsw];
  for (int i = 0; i < nsd; ++i) left_edges_rd[i] = 1e-6 * pow(10, -3 + i * .1 );  //TODO * si::metres
  for (int i = 0; i < nsw; ++i) left_edges_rw[i] = 1e-6 *  5 * pow(10,  0 + i * .025);  //TODO * si::metres
 
  notice_macro("opening netCDF file")
  NcFile nf("tmp/out.nc", NcFile::read);

  // be sure you're plotting for sdm simulation
  string micro, cmdline;
  nf.getAtt("command_line").getValues(cmdline); 
  if (cmdline.find("--eqs todo_sdm") != string::npos)
    micro = "sdm";
  else 
    assert(false);

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

//  for (size_t t = 0 ; t < nt; t + dt_out / si::seconds)
//  {
    Array<real_t, 2> tmp(nx, ny);
    for (int i = 0; i < nsw; ++i)
    {
      ostringstream name;
      name << "mw_0_" << i;
      nf.getVar(name.str()).getVar(start({focus_t,0,0,0}), count({1,nx,ny,1}), tmp.data());

      focus[left_edges_rw[i] / 1e-6] = sum(tmp(
	blitz::Range(focus_x-1, focus_x+1),
	blitz::Range(focus_y-1, focus_y+1)
      )) / 9 / 1e6;
    }

    string ext = "png";
    notice_macro("setting-up plot parameters");
    Gnuplot gp;

    gp << "reset" << endl;

    if (ext == "png")
      gp << "set term png enhanced size 1200," << 700 << endl;
    else if (ext == "eps")
      gp << "set term postscript size 40cm, 25cm solid enhanced color font 'Helvetica, 12'" << endl;
    else assert(false);

    gp << "set output 'tmp/focus_" << focus_t << ".png'" << endl;

  //  gp << "set title 'x/dx=" << focus_x << "+/-1    y/dy=" << focus_y << "+/-1     t/dt=" << t << endl;
    gp << "set logscale xy" << endl;
    gp << "set xrange [5:50]" << endl;
  //  gp << "set xtics (";
  //  for (int i = 0; i < left_edges.size()-1; i+=9)
  //    gp << (i!=0?",":"") << "\"" << format("%4.3g") % real_t(left_edges[i] / si::metres / real_t(1e-6)) << "\" " << i;
  //  gp << ")" << endl;
    gp << "set xlabel 'particle (dry) radius [{/Symbol m}m]'" << endl;
    gp << "set yrange [.001:100]" << endl;
    gp << "set ylabel 'particle concentration (per bin) [1/cm^3]'" << endl;
    gp << "set grid" << endl;
    gp << "plot '-' with steps lw 3 notitle" << endl;
    gp.send(focus);
  //}
}

