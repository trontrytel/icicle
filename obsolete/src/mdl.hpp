/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#pragma once
#include "opt_grd.hpp"
#include "opt_adv.hpp"
#include "opt_vel.hpp"
#include "opt_slv.hpp"
#include "opt_out.hpp"
#include "opt_ini.hpp"
#include "opt_eqs.hpp"
#include "opt_stp.hpp"

template <typename real_t>
void mdl(const po::variables_map &vm, const string &cmdline) 
{
  // grid choice
  notice_macro("parsing options: grd...")
  unique_ptr<grd<real_t>> grid(opt_grd<real_t>(vm));

  // advection scheme choice
  notice_macro("parsing options: adv...")
  unique_ptr<adv<real_t>> advsch(opt_adv<real_t>(vm));

  // velocity field choice
  notice_macro("parsing options: vel...")
  unique_ptr<vel<real_t>> velocity(opt_vel<real_t>(vm, *grid));

  // initial condition
  notice_macro("parsing options: ini...")
  unique_ptr<ini<real_t>> intcond(opt_ini<real_t>(vm, *grid));

  // equations
  notice_macro("parsing options: eqs...")
  unique_ptr<eqs<real_t>> eqsys(opt_eqs<real_t>(vm, *grid, *intcond, *velocity));

  // grouping all above into a single set-up object
  notice_macro("parsing options: stp...")
  unique_ptr<stp<real_t>> setup(opt_stp<real_t>(vm, 
    *advsch, 
    *velocity, 
    *intcond, 
    *grid, 
    *eqsys
  ));

  // output choice
  notice_macro("parsing options: out...")
  unique_ptr<out<real_t>> output(opt_out<real_t>(vm, *setup, cmdline));

  // solver choice
  notice_macro("parsing options: slv...")
  unique_ptr<slv<real_t>> solver(opt_slv(vm, *setup, *output));

  // integration
  notice_macro("parsing options: done.")
  solver->integ_loop(); 
}
