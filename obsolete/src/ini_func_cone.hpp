/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @brief contains definition of the ini_func_cone class
 */
#pragma once
#include "ini_func.hpp" 

/// @brief a cone-shaped 2D scalar field from the rotating cone test (Smolarkiewicz et al. 1983)
template <typename real_t>
class ini_func_cone : public ini_func<real_t>
{
  private: const quantity<si::length, real_t> h, x0, z0, r, h0;

  public: ini_func_cone(
    const grd<real_t> &grid,
    const quantity<si::length, real_t> &h,
    const quantity<si::length, real_t> &x0,
    const quantity<si::length, real_t> &z0,
    const quantity<si::length, real_t> &r,
    const quantity<si::length, real_t> &h0
  );

  public: quantity<si::dimensionless, real_t> psi(
    const quantity<si::length, real_t> &x,
    const quantity<si::length, real_t> &,
    const quantity<si::length, real_t> &z
  ) const;
};
