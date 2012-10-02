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


int main()
{
  string opts = wkc::icmw8_case1::create_files<real_t>(wkc::icmw8_case1::bulk);
 
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
