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

template <typename real_t>
class stp : root
{
  // that should probably be the only place where public fields are used
  // TODO: but maybe better offer public const-ref-returning methods?
  public: adv<real_t> *fllbck, *advsch;
  public: vel<real_t> *velocity;
  public: ini<real_t> *intcond;
  public: grd<real_t> *grid;

  public: int nx, ny, nz, nout;
  public: unsigned long nt;
  public: quantity<si::time, real_t> dt;

  public: stp(
    adv<real_t> *fllbck, 
    adv<real_t> *advsch, 
    vel<real_t> *velocity,
    ini<real_t> *intcond,
    grd<real_t> *grid,
    int nx, int ny, int nz, 
    quantity<si::time, real_t> dt_out, 
    quantity<si::time, real_t> t_max,
    quantity<si::time, real_t> dt_arg // if zero then calculate it automagically
  ) 
    : fllbck(fllbck), advsch(advsch), 
      velocity(velocity), 
      intcond(intcond), 
      grid(grid),
      nx(nx), ny(ny), nz(nz)
  { 
    bool auto_dt = (dt_arg/si::seconds == 0);
    dt = auto_dt ? real_t(1) * si::seconds : dt_arg;

    // TODO! merge with the logic from slv ctor and move into grid
    int halo = (advsch->stencil_extent() -1) / 2;
    mtx::rng
      ix = grid->rng_vctr(0, nx - 1, halo),
      jx = grid->rng_sclr(0, ny - 1, halo),
      kx = grid->rng_sclr(0, nz - 1, halo);
    mtx::rng
      iy = grid->rng_sclr(0, nx - 1, halo),
      jy = grid->rng_vctr(0, ny - 1, halo),
      ky = grid->rng_sclr(0, nz - 1, halo);
    mtx::rng
      iz = grid->rng_sclr(0, nx - 1, halo),
      jz = grid->rng_sclr(0, ny - 1, halo),
      kz = grid->rng_vctr(0, nz - 1, halo);
    mtx::arr<real_t>
      Cx(ix, jx, kx),
      Cy(iy, jy, ky),
      Cz(iz, jz, kz);

    grid->populate_courant_fields(
      ix, jx, kx,
      iy, jy, ky,
      iz, jz, kz,
      &Cx, &Cy, &Cz,
      velocity, dt,
      nx, ny, nz
    );
    
    real_t cmax = max(sqrt(pow2(Cx) + pow2(Cy) + pow2(Cz))); // TODO: check if that's a correct way to calculate it?

    if (auto_dt)
    {
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
    }
    else
    {
      nout = dt_out / dt;
      if (real_t(nout) * dt != dt_out) 
        error_macro("dt_out is not an integer multiple of dt");
      if (cmax > advsch->courant_max())
        warning_macro("chosen dt results in too big Courant number: max(C) = " << cmax << " > " << advsch->courant_max());
      if (cmax < advsch->courant_min())
        //donor cell schemes work better for big courant numbers
        warning_macro("chosen dt results in too small Courant number: min(C) = " << cmax << " < " << advsch->courant_min());
    }
    nt = t_max / dt;
  }
};

#endif
