/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date January 2012 - February 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 @  @brief contains definition of the eqs_shallow_water class - a system of 2D shallow-water equations
 */
#pragma once

#include "eqs.hpp"
#include "rhs_explicit.hpp"
#include "grd.hpp"
#include "phc/phc.hpp"

template <typename real_t>
class eqs_shallow_water : public eqs<real_t> 
{
  // TODO: these should become obsolete when psi becomes a map!
  private: struct params
  {
    quantity<si::length, real_t> h_unit;
    quantity<multiply_typeof_helper<si::velocity, si::length>::type, real_t> q_unit;
    quantity<si::dimensionless, real_t> dHdxy_unit;
    int idx_h;
  };  
  private: params par;

  private: template <int di, int dj> class forcings;

  // ctor
  public: eqs_shallow_water(const grd<real_t> &grid);
};
