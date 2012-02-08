/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date December 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 *  @section DESCRIPTION
 *    A simple container for storing simulation set-up elements, i.e.
 *    the advection scheme, the velocity field etc
 */
#ifndef STP_HPP
#  define STP_HPP

#  include "adv.hpp"
#  include "out.hpp"
#  include "vel.hpp"
#  include "ini.hpp"
#  include "grd.hpp"
#  include "eqs.hpp"

template <typename real_t>
struct stp : root
{
  adv<real_t> *fllbck, *advsch;
  vel<real_t> *velocity;
  ini<real_t> *intcond;
  grd<real_t> *grid;
  eqs<real_t> *equations;

  int nout;
  unsigned long nt;
  quantity<si::time, real_t> dt;

  stp(
    adv<real_t> *fllbck, 
    adv<real_t> *advsch, 
    vel<real_t> *velocity,
    ini<real_t> *intcond,
    grd<real_t> *grid,
    eqs<real_t> *equations,
    quantity<si::time, real_t> dt_out, 
    quantity<si::time, real_t> t_max
  ) 
    : fllbck(fllbck), advsch(advsch), 
      velocity(velocity), 
      intcond(intcond), 
      grid(grid),
      equations(equations)
  { 
    int halo = (advsch->stencil_extent() -1) / 2;
    mtx::idx_ijk ijk(
      mtx::rng(0, grid->nx() - 1),
      mtx::rng(0, grid->ny() - 1),
      mtx::rng(0, grid->nz() - 1)
    );
    mtx::arr<real_t>
      Cx(grid->rng_vctr_x(ijk, halo)),
      Cy(grid->rng_vctr_y(ijk, halo)),
      Cz(grid->rng_vctr_z(ijk, halo));
    grid->populate_courant_fields(&Cx, &Cy, &Cz, velocity, dt);
    
    real_t cmax = max(sqrt(pow2(Cx) + pow2(Cy) + pow2(Cz))); // TODO: check if that's a correct way to calculate it?

    dt = dt_out;
    nout = 1;
    while (cmax * dt / si::seconds > advsch->courant_max()) //TODO ? some limit to this loop
    {  
      dt = dt_out / real_t(++nout);
    }
    if (cmax*dt/si::seconds < advsch->courant_min()) 
      error_macro("failed to calculate a reasonable time step") 
      //TODO print some suggestions for t_out
    if (!velocity->is_constant())  
      warning_macro ("velocity field is not const -> calculated time step may change") 
    nt = t_max / dt;
  }

  stp(
    adv<real_t> *fllbck, 
    adv<real_t> *advsch, 
    vel<real_t> *velocity,
    ini<real_t> *intcond,
    grd<real_t> *grid,
    eqs<real_t> *equations,
    int nout, 
    quantity<si::time, real_t> dt, 
    unsigned long nt
  ) 
    : fllbck(fllbck), advsch(advsch), 
      velocity(velocity), 
      intcond(intcond), 
      grid(grid),
      equations(equations),
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
