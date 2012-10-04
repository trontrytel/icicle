/** @file
 *  @example wkc_icmw8_case1_t0/wkc_icmw8_case1_t0_test.cpp
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @author Sylwester Arabas <sarabas@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date October 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @section DESCRIPTION
 *    Tests the initial values for the ICMW case 1 set-up
 */

#include <map>
using std::map;

#include <boost/filesystem.hpp>

#include <boost/units/systems/si.hpp>
#include <boost/units/io.hpp>
namespace si = boost::units::si;
using boost::units::quantity;
using boost::units::divide_typeof_helper;
using boost::units::detail::get_value;

#define GNUPLOT_ENABLE_BLITZ
#include <gnuplot-iostream/gnuplot-iostream.h>
using std::endl;

#include "../../src/wkc/wkc_icmw8_case1.hpp"

typedef vector<size_t> start;
typedef vector<size_t> count;

typedef float real_t; // TODO: chose a supported one
using wkc::icmw8_case1::micro_str;
using wkc::icmw8_case1::micro_t;
using wkc::icmw8_case1::bulk;
using wkc::icmw8_case1::sdm;
using wkc::icmw8_case1::mm;

int main()
{
  list<enum micro_t> micros({bulk, sdm, mm});

  // initialising gnuplot
  Gnuplot gp;
  gp << "set term svg enhanced" << endl;
  gp << "set output 'tmp.svg'" << endl; 
  gp << "set xlabel 'r_l [g/kg]'" << endl;
  gp << "set key bottom right" << endl;
  gp << "set grid" << endl;
  gp << "plot 0 notitle";
  for (micro_t &micro : micros)
    gp << ",'-' using ($1*1000):0 with linespoints t \"r_l (" << micro_str[micro] << ")\"";
  gp << endl;

  // loop over types of microphysics
  for (micro_t &micro : micros)
  {
    // creating the input files
    boost::filesystem::remove_all("tmp_" + micro_str[micro]);
    boost::filesystem::create_directory("tmp_" + micro_str[micro]);
    string opts = wkc::icmw8_case1::create_files<real_t>(micro, 
      "tmp_" + micro_str[micro] + "/rho.nc", 
      "tmp_" + micro_str[micro] + "/ini.nc"
    );

    // checking ini files
    NcFile nfini("tmp_" + micro_str[micro] + "/ini.nc", NcFile::read);
    int ny = nfini.getDim("Y").getSize();
    Array<real_t, 1> rhod(ny);
    nfini.getVar("rhod").getVar(rhod.data());

    // sanity check
    if (min(rhod) < .9 || max(rhod) > 1.3) error_macro("the air seems too dense...");

    // TODO: check e.g. if really hydrostatic

    // running the simulations
    ostringstream cmd;
    cmd << "../../icicle " 
      << opts
      << " --adv mpdata"
      << " --dt_out 6"
      << " --t_max 6"
      << " --out netcdf"
      << " --out.netcdf.file tmp_" << micro_str[micro] << "/out.nc"
      << " --slv serial";

    if (micro == bulk) cmd 
      << " --eqs.todo_bulk.cevp true"
      << " --eqs.todo_bulk.conv false"
      << " --eqs.todo_bulk.clct false"
      << " --eqs.todo_bulk.sedi false"
      << " --eqs.todo_bulk.revp false";
    else if (micro == mm) cmd 
      << " --eqs.todo_mm.act   true"
      << " --eqs.todo_mm.cond  true"
      << " --eqs.todo_mm.autoc false"
      << " --eqs.todo_mm.acc   false"
      << " --eqs.todo_mm.self  false"
      << " --eqs.todo_mm.turb  false"
      << " --eqs.todo_mm.sedi  false";
    else if (micro == sdm) cmd 
      << " --eqs.todo_sdm.adve true"
      << " --eqs.todo_sdm.cond true"
      << " --eqs.todo_sdm.coal false"
      << " --eqs.todo_sdm.sedi false"
      << " --eqs.todo_sdm.chem false";

    if (EXIT_SUCCESS != system(cmd.str().c_str()))
      error_macro("model run failed: " << cmd.str())

    {
      NcFile nfout("tmp_" + micro_str[micro] + "/out.nc", NcFile::read);
      int nx = nfout.getDim("X").getSize();
      if (ny != nfout.getDim("Y").getSize()) error_macro("ini:ny != out:ny");
      Array<real_t, 2> rhod_rl(nx, ny);

      int t = 1; // letting the non-bulk models to deal with supersaturation
      
      if (micro == sdm)
      {
        nfout.getVar("m_3").getVar(start({t,0,0,0}), count({1,nx,ny,1}), rhod_rl.data());
        rhod_rl *= real_t(4./3.) * phc::pi<real_t>() * phc::rho_w<real_t>() * si::cubic_metres / si::kilograms;
      }
      else
      {
        nfout.getVar("rhod_rl").getVar(start({t,0,0,0}), count({1,nx,ny,1}), rhod_rl.data());
      }

      Array<real_t, 1> rl_mean(ny);
      rl_mean = 0;
      for (int y = 0; y < ny; ++y)
      {
        for (int x = 0; x < nx; ++x) 
          rl_mean(y) += rhod_rl(x, y);
        rl_mean(y) /= nx * rhod(y);
      }
      gp.send(rl_mean);
    }
  }
}
