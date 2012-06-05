/* @file
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date March 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#include <vector>
using std::vector;

#include <map>
using std::map;

#include <iostream>
using std::ostringstream;
using std::cerr;
using std::endl;
#define notice_macro(msg) { cerr << msg << endl; }

#include <boost/units/systems/si.hpp>
#include <boost/units/io.hpp>
namespace si = boost::units::si;
using boost::units::quantity;
using boost::units::divide_typeof_helper;
using boost::units::detail::get_value;

#define GNUPLOT_ENABLE_BLITZ
#include <gnuplot-iostream/gnuplot-iostream.h>

#include "../../src/cmn.hpp"
#include "../../src/phc_theta.hpp"
#include "../../src/phc_terminal_vel.hpp"

typedef float real_t;

int main()
{
  quantity<si::temperature, real_t> T=real_t(293.)*si::kelvin;
  quantity<si::mass_density, real_t> rhoa=real_t(1.23)*si::kilograms/si::cubic_metres;
  quantity<si::length, real_t> r=real_t(1.*1e-6)*si::metres;

  quantity<si::length, real_t> radius;
  quantity<si::velocity, real_t> term_vel;
 
  cout << 1e2 << endl;
  cout << 10e2 << endl;

  map< quantity<si::length, real_t>, quantity<si::velocity, real_t> > plot_data;
  for(int i=0; i<=45; i++){
    radius = r*real_t(i*100) ; 
//    cout << radius << endl;
    term_vel = phc::vt(r * real_t(i * 100),T, rhoa) ;
//    cout << term_vel << endl;
    plot_data[radius * real_t(2*1e6)] = term_vel * real_t(100);
  }
 
  Gnuplot gp;
  gp << "set xlabel 'drop diameter [um]'" << endl;
  gp << "set ylabel 'terminal velocity [cm/s]'" << endl;
  gp << "set term postscript" << endl;
  gp << "set output 'term_vel_test.ps'" << endl; 
  gp << "plot '-' u 1:3" << endl;
  gp.send(plot_data);

}
