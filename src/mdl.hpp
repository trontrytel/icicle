/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef MODEL_HPP
#  define MODEL_HPP

#  include "opt_grd.hpp"
#  include "opt_adv.hpp"
#  include "opt_vel.hpp"
#  include "opt_slv.hpp"
#  include "opt_out.hpp" // has to be included after slv/MPI
#  include "opt_ini.hpp"

template <typename real_t>
void mdl(const po::variables_map &vm, const string &options) 
{
  // some key parameters (TODO: move from here!)
  if (
    !vm.count("nx") || !vm.count("ny") || !vm.count("nz") || 
    !vm.count("nt") || !vm.count("dt")
  )
    error_macro("nx, ny, nz, nt and dt options are mandatory")
  quantity<si::time, real_t> 
    dt = boost::lexical_cast<real_t>(vm["dt"].as<string>()) * si::seconds;
  int 
    nx = vm["nx"].as<int>(),
    ny = vm["ny"].as<int>(),
    nz = vm["nz"].as<int>();
  unsigned long
    nt = vm["nt"].as<unsigned long>();

  // sanity checks
  if (nx <= 0 || ny <= 0 || nz <= 0) error_macro("nx, ny, nz must all be >= 0") 

  // grid choice
  auto_ptr<grd<real_t> > grid(opt_grd<real_t>(vm));

  // advection scheme choice
  auto_ptr<adv<real_t> > advsch, fllbck;
  {
    adv<real_t> *advschp, *fllbckp;
    opt_adv<real_t>(vm, &fllbckp, &advschp, grid.get());
    advsch.reset(advschp);
    if (fllbckp != NULL) fllbck.reset(fllbckp);
  }

  // output choice
  auto_ptr<out<real_t> > output(opt_out<real_t>(vm, grid.get(), nx, ny, nz, options));

  // velocity choice
  auto_ptr<vel<real_t> > velocity(opt_vel<real_t>(vm, grid.get(), nx, ny, nz));

  // initial condition
  auto_ptr<ini<real_t> > intcond(opt_ini<real_t>(vm, grid.get()));

  // grouping all above into a single set-up object
  auto_ptr<stp<real_t> > setup(new stp<real_t>(
    fllbck.get(), advsch.get(), 
    output.get(), 
    velocity.get(),
    intcond.get(),
    grid.get()
  ));

  // solver choice
  auto_ptr<slv<real_t> > solver(opt_slv(vm, setup.get(), nx, ny, nz, dt));

  // integration
  solver->integ_loop(nt, dt);
}
#endif
