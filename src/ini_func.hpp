/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date January 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @brief contains definition of the ini_func class - a base class for functional-form specified initial conditions
 */
#pragma once
#include "ini.hpp"
#include "grd.hpp"

/// @brief a base class for functional-form specified initial conditions (still pure virtual, applicable to single-equation systems only)
template <typename real_t>
class ini_func : public ini<real_t>
{
  private: const grd<real_t> *grid;
  public: ini_func(const grd<real_t> &grid)
    : grid(&grid)
  { }

  private: virtual quantity<si::dimensionless, real_t> psi(
    const quantity<si::length, real_t> &x,
    const quantity<si::length, real_t> &y,
    const quantity<si::length, real_t> &z
  ) const = 0;

  public: virtual void populate_scalar_field(
    const string &varname,
    const mtx::idx &ijk,
    mtx::arr<real_t> &data
  ) const
  {
    if (varname != "psi") 
      error_macro("specification of init. cond. with a function works for single-equation systems only")

    for (int i = ijk.lbound(mtx::i); i <= ijk.ubound(mtx::i); ++i)
      for (int j = ijk.lbound(mtx::j); j <= ijk.ubound(mtx::j); ++j)
        for (int k = ijk.lbound(mtx::k); k <= ijk.ubound(mtx::k); ++k)
          data(i,j,k) = psi(grid->x(i,j,k), grid->y(i,j,k), grid->z(i,j,k));
  }
};
