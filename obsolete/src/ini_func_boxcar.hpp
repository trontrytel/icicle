/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011 - January 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @brief contains definition of the ini_func_boxcar class
 */
#pragma once
#include "ini_func.hpp" 

/// @brief a multidimensional boxcar shape definition 
template <typename real_t>
class ini_func_boxcar : public ini_func<real_t>
{
  private: const quantity<si::length, real_t> ax, bx, ay, by, az, bz;
  private: const quantity<si::dimensionless, real_t> A, A0;

  public: ini_func_boxcar(
    const grd<real_t> &grid,
    const quantity<si::length, real_t> &ax,
    const quantity<si::length, real_t> &bx,
    const quantity<si::length, real_t> &ay,
    const quantity<si::length, real_t> &by,
    const quantity<si::length, real_t> &az,
    const quantity<si::length, real_t> &bz,
    const quantity<si::dimensionless, real_t> &A,
    const quantity<si::dimensionless, real_t> &A0
  );

  public: quantity<si::dimensionless, real_t> psi(
    const quantity<si::length, real_t> &x,
    const quantity<si::length, real_t> &y,
    const quantity<si::length, real_t> &z
  ) const;
};
