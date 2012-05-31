/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef MODEL_HPP
#  define MODEL_HPP

#  include "opt_grd.hpp"
#  include "opt_adv.hpp"
#  include "opt_vel.hpp"
#  include "opt_slv.hpp"
#  include "opt_out.hpp"
#  include "opt_ini.hpp"
#  include "opt_eqs.hpp"
#  include "opt_stp.hpp"

template <typename real_t>
void mdl(const po::variables_map &vm, const string &cmdline) 
{
  // grid choice
  cerr << "-- init: parsing options: grd..." << endl;
  unique_ptr<grd<real_t>> grid(opt_grd<real_t>(vm));

  // advection scheme choice
  cerr << "-- init: parsing options: adv..." << endl;
  unique_ptr<adv<real_t>> advsch(opt_adv<real_t>(vm, grid.get()));

  // velocity field choice
  cerr << "-- init: parsing options: vel..." << endl;
  unique_ptr<vel<real_t>> velocity(opt_vel<real_t>(vm, *grid));

  // initial condition
  cerr << "-- init: parsing options: ini..." << endl;
  unique_ptr<ini<real_t>> intcond(opt_ini<real_t>(vm, *grid));

  // equations
  cerr << "-- init: parsing options: eqs..." << endl;
  unique_ptr<eqs<real_t>> eqsys(opt_eqs<real_t>(vm, *grid, *intcond, *velocity));

  // grouping all above into a single set-up object
  cerr << "-- init: parsing options: stp..." << endl;
  unique_ptr<stp<real_t>> setup(opt_stp<real_t>(vm, 
    advsch.get(), 
    velocity.get(), 
    intcond.get(), 
    grid.get(), 
    eqsys.get()
  ));

  // output choice
  cerr << "-- init: parsing options: out..." << endl;
  unique_ptr<out<real_t>> output(opt_out<real_t>(vm, setup.get(), cmdline));

  // solver choice
  cerr << "-- init: parsing options: slv..." << endl;
  unique_ptr<slv<real_t>> solver(opt_slv(vm, setup.get(), output.get()));

  // integration
  cerr << "-- init: parsing options: done." << endl;
  solver->integ_loop(); 
}
#endif
