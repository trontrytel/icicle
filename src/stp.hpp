/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date December 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @brief contains definition of the @ref stp class - a container for simulation set-up elements
 */
#ifndef STP_HPP
#  define STP_HPP

#  include "adv.hpp"
#  include "out.hpp"
#  include "vel.hpp"
#  include "ini.hpp"
#  include "grd.hpp"
#  include "eqs.hpp"

/// @brief A simple container for storing simulation set-up elements, i.e.
///        advection scheme, velocity field, initial condition, grid, equation system, etc
template <typename real_t>
struct stp 
{
  const adv<real_t> &advsch;
  const vel<real_t> &velocity;
  const ini<real_t> &intcond;
  const grd<real_t> &grid;
  eqs<real_t> &eqsys; // cannot be const as the super-droplet state is kept herein (TODO?)

  int nout;
  unsigned long nt;
  quantity<si::time, real_t> dt;

  stp(
    const adv<real_t> &advsch, 
    const vel<real_t> &velocity,
    const ini<real_t> &intcond,
    const grd<real_t> &grid,
    eqs<real_t> &eqsys,
    const quantity<si::time, real_t> dt_out, 
    const quantity<si::time, real_t> t_max
  ) 
    : advsch(advsch), 
      velocity(velocity), 
      intcond(intcond), 
      grid(grid),
      eqsys(eqsys)
  { 
    // TODO: it does not work for non-constant velocities as of now!
    int halo = (advsch.stencil_extent() -1) / 2;
    mtx::idx_ijk ijk(
      mtx::rng(0, grid.nx() - 1),
      mtx::rng(0, grid.ny() - 1),
      mtx::rng(0, grid.nz() - 1)
    );
    mtx::arr<real_t>
      Cx(grid.rng_vctr_x(ijk, halo)),
      Cy(grid.rng_vctr_y(ijk, halo)),
      Cz(grid.rng_vctr_z(ijk, halo));

    dt = dt_out;
    nout = 1;

    velocity.populate_courant_fields(-1, -1, &Cx, &Cy, &Cz, dt, NULL, NULL, NULL); // TODO: only if vel->is_constant() !!!
    // Cx, Cy and Cz dimensions are not the same with Arakawa-C grid!
    //real_t cmax = max(sqrt(pow2(Cx) + pow2(Cy) + pow2(Cz))); // TODO: check if that's a correct way to calculate it?
    real_t cmax = max(abs(Cx)) + max(abs(Cy)) + max(abs(Cz));
//cerr << max(abs(Cx)) << " + " << max(abs(Cy)) << " + " << max(abs(Cz)) << endl;

    while (cmax * dt / si::seconds > advsch.courant_max()) //TODO ? some limit to this loop
    {  
      dt = dt_out / real_t(++nout);
    }
    if (cmax*dt/si::seconds < advsch.courant_min()) 
      error_macro("failed to calculate a reasonable time step for" << endl
        << "t_max=" << t_max << endl
        << "dt_out=" << dt_out << endl
        << "cmax = " << cmax
      ) 
    if (!velocity.is_constant())  
      warning_macro("velocity field is not const -> calculated time step may change") 
    nt = t_max / dt;
  }

  stp(
    const adv<real_t> &advsch, 
    const vel<real_t> &velocity,
    const ini<real_t> &intcond,
    const grd<real_t> &grid,
    eqs<real_t> &eqsys,
    const int nout, 
    const quantity<si::time, real_t> dt, 
    const unsigned long nt
  ) 
    : advsch(advsch), 
      velocity(velocity), 
      intcond(intcond), 
      grid(grid),
      eqsys(eqsys),
      nt(nt), dt(dt), nout(nout)
  {
  /*       if (cmax > advsch->courant_max())  TODO 
        warning_macro("chosen dt results in too big Courant number: max(C) = " << cmax << " > " << advsch->courant_max());
      if (cmax < advsch->courant_min())
        //donor cell schemes work better for big courant numbers
        warning_macro("chosen dt results in too small Courant number: min(C) = " << cmax << " < " << advsch->courant_min());
  */
  }  

};

#endif
