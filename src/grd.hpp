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
    Array<quantity<si::dimensionless, real_t>, 3> *Cx, 
    Array<quantity<si::dimensionless, real_t>, 3> *Cy, 
    Array<quantity<si::dimensionless, real_t>, 3> *Cz, 
    vel<real_t> *velocity,
    quantity<si::time, real_t> dt
  ) 
  {
    for (int i = ir.first(); i <= ir.last(); ++i)
      for (int j = jr.first(); j <= jr.last(); ++j)
        for (int k = kr.first(); k <= kr.last(); ++k)
        {
          if (true)
            (*Cx)(i, j, k) = dt / dx() *
              velocity->u(u_x(i, j, k), (jr.first() != jr.last()) ? u_y(i, j, k) : 0, u_z(i, j, k));
          if (jr.first() != jr.last()) 
            (*Cy)(i, j, k) = dt / dy() *
              velocity->v(v_x(i, j, k), v_y(i, j, k), v_z(i, j, k));
          if (kr.first() != kr.last()) 
            (*Cz)(i, j, k) = dt / dz() *
              velocity->w(w_x(i, j, k), (jr.first() != jr.last()) ? w_y(i, j, k) : 0, w_z(i, j, k));
        }
  }
  
};
#endif
