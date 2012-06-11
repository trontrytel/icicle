/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date January 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#pragma once
#include "ini_func.hpp" 

/// @brief a multi-dimensional gaussian shape definition
template <typename real_t>
class ini_func_gauss : public ini_func<real_t>
{
  private: const quantity<si::dimensionless, real_t> A, A0;
  private: const quantity<si::length, real_t> x0, y0, z0, sx, sy, sz;

  public: ini_func_gauss(
    const grd<real_t> &grid,
    const quantity<si::dimensionless, real_t> &A,
    const quantity<si::dimensionless, real_t> &A0,
    const quantity<si::length, real_t> &x0,
    const quantity<si::length, real_t> &y0,
    const quantity<si::length, real_t> &z0,
    const quantity<si::length, real_t> &sx,
    const quantity<si::length, real_t> &sy,
    const quantity<si::length, real_t> &sz
  );

  public: quantity<si::dimensionless, real_t> psi(
    const quantity<si::length, real_t> &x,
    const quantity<si::length, real_t> &y,
    const quantity<si::length, real_t> &z
  ) const;
};
