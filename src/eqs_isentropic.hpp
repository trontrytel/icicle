/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date February 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 @  @brief contains definition of the eqs_isentropic class - a system of 3D isentropic equations
 */
#pragma once

#include "eqs.hpp"
#include "rhs_implicit.hpp"
#include "grd.hpp"
#include "phc/phc_theta.hpp"

/// @brief the 3D isentropic equations system
template <typename real_t>
class eqs_isentropic : public eqs<real_t> 
{
  // TODO: should become obsolete when psi becomes a map
  private: struct params
  {
    // TODO: add const qualifiers and initiliase here if possible
    quantity<si::pressure, real_t> p_unit, p_top;
    quantity<multiply_typeof_helper<si::velocity, si::pressure>::type, real_t> q_unit;
    quantity<si::temperature, real_t> theta_unit, theta_frst;
    quantity<multiply_typeof_helper<si::acceleration, si::length>::type, real_t> M_unit; 
    quantity<si::dimensionless, real_t> dHdxy_unit; 
    vector<int> idx_dp;
  };

  // nested classes
  private: class rayleigh_damping;
  private: template <int di, int dj> class montgomery_grad;

  // private fields
  private: params par;

  // ctor
  public: eqs_isentropic(const grd<real_t> &grid,
    int nlev, 
    quantity<si::pressure, real_t> p_top,
    quantity<si::temperature, real_t> theta_frst,
    int abslev, real_t absamp
  );
};
