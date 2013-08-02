/** @file
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date August 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#pragma once
#include "phc.hpp"

namespace phc
{
  namespace kelvin
  {
    // water surface tension (Petters and Kreidenweis 2007 //TODO - check)
    phc_declare_const_macro(sg_surf, 0.072, si::newton/si::metres);


    // Kelvin curvature parameter (see eq. 7 in Kvorostyanov and Curry 2006)
    phc_declare_funct_macro quantity<si::length, real_t> A(
      quantity<si::temperature, real_t> T
    )
    {
      return real_t(2.) * sg_surf<real_t>() / R_v<real_t>() / T / rho_w<real_t>();
    }

    phc_declare_funct_macro quantity<si::dimensionless, real_t> klvntrm(
      quantity<si::length, real_t> r,
      quantity<si::temperature, real_t> T
    )
    {
      return exp(real_t(2.) * sg_surf<real_t>() / r / R_v<real_t>() / rho_w<real_t>() / T);
    }
  }
};
