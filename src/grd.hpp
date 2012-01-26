/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef GRD_HPP
#  define GRD_HPP

#  include "cmn.hpp" // root class, error reporting
#  include "vel.hpp"
#  include "ini.hpp"
#  include "mtx.hpp"

template<typename real_t>
class grd : root
{
  public: virtual quantity<si::length, real_t> dx() = 0;
  public: virtual quantity<si::length, real_t> dy() = 0;
  public: virtual quantity<si::length, real_t> dz() = 0;

  public: virtual int nx() = 0;
  public: virtual int ny() = 0;
  public: virtual int nz() = 0;

  // returns ranges to be passed as constructors to arr
  // first and last reflect scalar indices
  public: virtual mtx::rng rng_sclr(int first, int last, int halo) = 0;
  public: virtual mtx::rng rng_vctr(int first, int last, int halo) = 0;

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
  public: virtual void populate_courant_fields(
    mtx::rng &ix, mtx::rng &jx, mtx::rng &kx,
    mtx::rng &iy, mtx::rng &jy, mtx::rng &ky,
    mtx::rng &iz, mtx::rng &jz, mtx::rng &kz,
    mtx::arr<real_t> *Cx, 
    mtx::arr<real_t> *Cy, 
    mtx::arr<real_t> *Cz, 
    vel<real_t> *velocity,
    quantity<si::time, real_t> dt,
    int nx, int ny, int nz
  ) 
  {
    quantity<si::length, real_t> 
      mx = real_t(nx) * dx(),
      my = real_t(ny) * dy(),
      mz = real_t(nz) * dz();
    for (int i = ix.first(); i <= ix.last(); ++i)
      for (int j = jx.first(); j <= jx.last(); ++j)
        for (int k = kx.first(); k <= kx.last(); ++k)
          (*Cx)(i, j, k) = dt / dx() *
            velocity->u(fmod(u_x(i, j, k), mx), fmod(u_y(i, j, k), my), fmod(u_z(i, j, k), mz));
    for (int i = iy.first(); i <= iy.last(); ++i)
      for (int j = jy.first(); j <= jy.last(); ++j)
        for (int k = ky.first(); k <= ky.last(); ++k)
          (*Cy)(i, j, k) = dt / dy() *
            velocity->v(fmod(v_x(i, j, k), mx), fmod(v_y(i, j, k), my), fmod(v_z(i, j, k), mz));
    for (int i = iz.first(); i <= iz.last(); ++i)
      for (int j = jz.first(); j <= jz.last(); ++j)
        for (int k = kz.first(); k <= kz.last(); ++k)
          (*Cz)(i, j, k) = dt / dz() *
            velocity->w(fmod(w_x(i, j, k), mx), fmod(w_y(i, j, k), my), fmod(w_z(i, j, k), mz));
  }

  public: virtual void populate_scalar_field(
    const mtx::rng &ii, const mtx::rng &jj, const mtx::rng &kk,
    mtx::arr<real_t> *psi, ini<real_t> *intcond
  )
  {
    for (int i = ii.first(); i <= ii.last(); ++i)
      for (int j = jj.first(); j <= jj.last(); ++j)
        for (int k = kk.first(); k <= kk.last(); ++k)
          (*psi)(i,j,k) = intcond->psi(x(i,j,k), y(i,j,k), z(i,j,k));
  }
};
#endif
