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

#include <list>
using std::list;

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

typedef float real_t; // TODO: chose a supported one
using wkc::icmw8_case1::micro_t;
using wkc::icmw8_case1::bulk;
using wkc::icmw8_case1::sdm;
using wkc::icmw8_case1::mm;

//TODO check ini files

int main()
{
  for (micro_t &micro : list<enum micro_t>({bulk, sdm, mm}))
  {
    map<enum micro_t, string> micro_str({{bulk, "bulk"}, {sdm, "sdm"}, {mm, "mm"}});

    boost::filesystem::remove_all("tmp_" + micro_str[micro]);
    boost::filesystem::create_directory("tmp_" + micro_str[micro]);
    string opts = wkc::icmw8_case1::create_files<real_t>(micro, 
      "tmp_" + micro_str[micro] + "/rho.nc", 
      "tmp_" + micro_str[micro] + "/ini.nc"
    );

    ostringstream cmd;
    cmd
      << "../../icicle " << opts
        << " --adv mpdata"
        << " --dt_out 3"
        << " --t_max 3"
        << " --out netcdf"
          << " --out.netcdf.file tmp_" << micro_str[micro] << "/out.nc"
       << " --slv serial"
       << " --eqs todo_" + micro_str[micro];
//     if (micro == sdm) 

    if (EXIT_SUCCESS != system(cmd.str().c_str()))
      error_macro("model run failed: " << cmd.str())
  }
/*
  map< quantity<si::length, real_t>, quantity<si::velocity, real_t> > plot_data;
  for(int i=0; i<=45; i++){
    radius = r*real_t(i*100) ; 
    term_vel = phc::vt(r * real_t(i * 100),T, rhoa) ;
    plot_data[radius * real_t(2*1e6)] = term_vel * real_t(100);
  }
 
  Gnuplot gp;
  gp << "set xlabel 'drop diameter [um]'" << endl;
  gp << "set ylabel 'terminal velocity [cm/s]'" << endl;
  gp << "set term svg" << endl;
  gp << "set output 'term_vel_test.svg'" << endl; 
  gp << "plot '-' u 1:3" << endl;
  gp.send(plot_data);
*/
}
