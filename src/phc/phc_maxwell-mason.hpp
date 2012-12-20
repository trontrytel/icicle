/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date March 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#pragma once
#include "phc.hpp"

namespace phc
{
  phc_declare_const_macro(D, 2.21e-5, si::square_metres / si::seconds) // water diffusivity in air
  phc_declare_const_macro(K, 2.4e-2, si::joules / si::metres / si::seconds / si::kelvins) // thermal conductivity

  phc_declare_funct_macro quantity<divide_typeof_helper<si::area, si::time>::type, real_t> rdrdt(
    const quantity<si::mass_density, real_t> rho_v, // ambient water vapour density
    const quantity<si::temperature, real_t> T, // ambient temperature
    const quantity<si::dimensionless, real_t> a_w, // water activity
    const quantity<si::dimensionless, real_t> klvntrm // the Kelvin term
  )
  {
    quantity<si::pressure, real_t> p_vs = phc::p_vs<real_t>(T);
    quantity<divide_typeof_helper<si::energy, si::mass>::type, real_t> l_v = phc::l_v<real_t>(T);
    quantity<divide_typeof_helper<si::energy, si::mass>::type, real_t> RvT = R_v<real_t>() * T;
    return 
      ( // vapour density difference (ambient - drop surface)
        rho_v - (p_vs / RvT * a_w) * klvntrm
      ) 
      / phc::rho_w<real_t>() // water density
      / (
        real_t(1) / D<real_t>() 
        +
        l_v / K<real_t>() / T * p_vs / RvT * (l_v / RvT - real_t(1))
      )
    ;
  }
};
