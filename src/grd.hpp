/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef GRD_HPP
#  define GRD_HPP

#  include "config.hpp" // USE_* defines
#  include "common.hpp" // root class, error reporting
#  include "vel.hpp"
#  include "ini.hpp"

template<typename real_t>
class grd : root
{
  public: virtual quantity<si::length, real_t> dx() = 0;
  public: virtual quantity<si::length, real_t> dy() = 0;
  public: virtual quantity<si::length, real_t> dz() = 0;

  // returns ranges to be passed as constructors to arr (e.g.)
  // first and last concern scalar indices
//  public: virtual Range rng_sclr(int first, int last) = 0;
  public: virtual Range rng_vctr(int first, int last) = 0;

  public: virtual quantity<si::length, real_t> x(int i, int j, int k) = 0;
  public: virtual quantity<si::length, real_t> y(int i, int j, int k) = 0;
  public: virtual quantity<si::length, real_t> z(int i, int j, int k) = 0;

  public: virtual quantity<si::length, real_t> u_x(int i, int j, int k) = 0;
  public: virtual quantity<si::length, real_t> u_y(int i, int j, int k) = 0;
  public: virtual quantity<si::length, real_t> u_z(int i, int j, int k) = 0;
  public: virtual quantity<si::length, real_t> v_x(int i, int j, int k) = 0;
  public: virtual quantity<si::length, real_t> v_y(int i, int j, int k) = 0;
  public: virtual quantity<si::length, real_t> v_z(int i, int j, int k) = 0;
  public: virtual quantity<si::length, real_t> w_x(int i, int j, int k) = 0;
  public: virtual quantity<si::length, real_t> w_y(int i, int j, int k) = 0;
  public: virtual quantity<si::length, real_t> w_z(int i, int j, int k) = 0;

  // ... Nothing is real and nothing to get hung about.
  // courant fields forever ...
  public: virtual void populate_courant_fields(Range &ir, Range &jr, Range &kr,
    Array<real_t, 3> *Cx, 
    Array<real_t, 3> *Cy, 
    Array<real_t, 3> *Cz, 
    vel<real_t> *velocity,
    quantity<si::time, real_t> dt
  ) 
  {
    for (int i = ir.first(); i <= ir.last(); ++i)
      for (int j = jr.first(); j <= jr.last(); ++j)
        for (int k = kr.first(); k <= kr.last(); ++k)
        {
          (*Cx)(i, j, k) = dt / dx() *
            velocity->u(u_x(i, j, k), u_y(i, j, k), u_z(i, j, k));
          (*Cy)(i, j, k) = dt / dy() *
            velocity->v(v_x(i, j, k), v_y(i, j, k), v_z(i, j, k));
          (*Cz)(i, j, k) = dt / dz() *
            velocity->w(w_x(i, j, k), w_y(i, j, k), w_z(i, j, k));
        }
  }

  public: virtual void populate_scalar_field(
    const Range &ii, const Range &jj, const Range &kk,
    Array<real_t, 3> *psi, ini<real_t> *intcond
  )
  {
    for (int i = ii.first(); i <= ii.last(); ++i)
      for (int j = jj.first(); j <= jj.last(); ++j)
        for (int k = kk.first(); k <= kk.last(); ++k)
          (*psi)(i,j,k) = intcond->psi(x(i,j,k), y(i,j,k), z(i,j,k));
  }
};
#endif
