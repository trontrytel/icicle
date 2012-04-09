/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011 - February 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef VEL_FUNC_HPP
#  define VEL_FUNC_HPP

#  include "vel.hpp" 

template <typename real_t>
class vel_func : public vel<real_t>
{
  public: bool is_constant() { return true; } 

  private: virtual quantity<si::velocity, real_t> u(
    const quantity<si::length, real_t> &x,
    const quantity<si::length, real_t> &y,
    const quantity<si::length, real_t> &z
  ) = 0;
  private: virtual quantity<si::velocity, real_t> v(
    const quantity<si::length, real_t> &x,
    const quantity<si::length, real_t> &y,
    const quantity<si::length, real_t> &z
  ) = 0;
  private: virtual quantity<si::velocity, real_t> w(
    const quantity<si::length, real_t> &x,
    const quantity<si::length, real_t> &y,
    const quantity<si::length, real_t> &z
  ) = 0;

  private: const grd<real_t> &grid;
  public: vel_func(const grd<real_t> &grid)
    : grid(grid)
  { }

  public: void populate_courant_fields(int,int,
    mtx::arr<real_t> *Cx, 
    mtx::arr<real_t> *Cy, 
    mtx::arr<real_t> *Cz, 
    quantity<si::time, real_t> dt,
    mtx::arr<real_t> *[],
    mtx::arr<real_t> *[],
    mtx::arr<real_t> *[]
  )   
  {
    quantity<si::length, real_t> 
      mx = real_t(grid.nx()) * grid.dx(),
      my = real_t(grid.ny()) * grid.dy(),
      mz = real_t(grid.nz()) * grid.dz();
    for (int i = Cx->lbound(mtx::i); i <= Cx->ubound(mtx::i); ++i)
      for (int j = Cx->lbound(mtx::j); j <= Cx->ubound(mtx::j); ++j)
        for (int k = Cx->lbound(mtx::k); k <= Cx->ubound(mtx::k); ++k)
          (*Cx)(i, j, k) = dt / grid.dx() *
            this->u(fmod(grid.u_x(i, j, k), mx), fmod(grid.u_y(i, j, k), my), fmod(grid.u_z(i, j, k), mz));
    for (int i = Cy->lbound(mtx::i); i <= Cy->ubound(mtx::i); ++i)
      for (int j = Cy->lbound(mtx::j); j <= Cy->ubound(mtx::j); ++j)
        for (int k = Cy->lbound(mtx::k); k <= Cy->ubound(mtx::k); ++k)
          (*Cy)(i, j, k) = dt / grid.dy() *
            this->v(fmod(grid.v_x(i, j, k), mx), fmod(grid.v_y(i, j, k), my), fmod(grid.v_z(i, j, k), mz));
    for (int i = Cz->lbound(mtx::i); i <= Cz->ubound(mtx::i); ++i)
      for (int j = Cz->lbound(mtx::j); j <= Cz->ubound(mtx::j); ++j)
        for (int k = Cz->lbound(mtx::k); k <= Cz->ubound(mtx::k); ++k)
          (*Cz)(i, j, k) = dt / grid.dz() *
            this->w(fmod(grid.w_x(i, j, k), mx), fmod(grid.w_y(i, j, k), my), fmod(grid.w_z(i, j, k), mz));
  }
};
#endif
