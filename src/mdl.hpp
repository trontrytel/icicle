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
  unique_ptr<grd<real_t> > grid(opt_grd<real_t>(vm));

  // advection scheme choice
  unique_ptr<adv<real_t> > advsch, fllbck;
  {
    adv<real_t> *advschp, *fllbckp;
    opt_adv<real_t>(vm, &fllbckp, &advschp, grid.get());
    advsch.reset(advschp);
    if (fllbckp != NULL) fllbck.reset(fllbckp);
  }

  // velocity field choice
  unique_ptr<vel<real_t> > velocity(opt_vel<real_t>(vm, grid.get()));

  // initial condition
  unique_ptr<ini<real_t> > intcond(opt_ini<real_t>(vm, *grid));

  // equations
  unique_ptr<eqs<real_t> > eqsys(opt_eqs<real_t>(vm, *grid));

  // grouping all above into a single set-up object
  unique_ptr<stp<real_t> > setup(opt_stp<real_t>(vm, 
    advsch.get(), 
    fllbck.get(), 
    velocity.get(), 
    intcond.get(), 
    grid.get(), 
    eqsys.get()
  ));

  // output choice
  unique_ptr<out<real_t> > output(opt_out<real_t>(vm, setup.get(), cmdline));

  // solver choice
  unique_ptr<slv<real_t> > solver(opt_slv(vm, setup.get(), output.get()));

  // integration
  solver->integ_loop(); 
}
#endif
