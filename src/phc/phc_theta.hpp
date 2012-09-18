/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date March 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @brief contains definition of thermodynamic relations related
 *    to the definition of potential temperature and the Exner function
 */
#pragma once
#include "phc.hpp"
#include "phc_const_cp.hpp" // TODO: option!

namespace phc
{
  // Exner function exponent for dry air
  phc_derived_const_macro(R_d_over_c_pd, R_d<real_t>() / c_pd<real_t>())

  // Exner function exponent for moist air
  phc_declare_funct_macro quantity<si::dimensionless, real_t> R_over_c_p(
    quantity<mixing_ratio, real_t> r
  )
  {
    return R<real_t>(r) / c_p<real_t>(r);
  }


  // Exner function for dry air
  phc_declare_funct_macro quantity<si::dimensionless, real_t> exner(
    quantity<si::pressure, real_t> p
  )
  {
    return pow(p / p_1000<real_t>(), R_d_over_c_pd<real_t>());
  }

  // Exner function for moist air
  phc_declare_funct_macro quantity<si::dimensionless, real_t> exner(
    quantity<si::pressure, real_t> p, 
    quantity<mixing_ratio, real_t> r
  )
  {
    return pow(p / p_1000<real_t>(), R_over_c_p<real_t>(r));
  }


  // temperature as a function theta, pressure and water vapour mixing ratio for moist air
  phc_declare_funct_macro quantity<si::temperature, real_t> T(
    quantity<si::temperature, real_t> th, // TODO: const
    quantity<si::pressure, real_t> p, // TODO: const
    quantity<mixing_ratio, real_t> r // TODO: const
  )
  {
    return th * exner<real_t>(p, r);
  }


  // pressure as a function of "theta times dry air density" and water vapour mixing ratio
  phc_declare_funct_macro quantity<si::pressure, real_t> p(
    const quantity<multiply_typeof_helper<si::mass_density, si::temperature>::type, real_t> rhod_th,
    quantity<mixing_ratio, real_t> r // TODO: const
  )
  {
    return p_1000<real_t>() * real_t(pow(
      (rhod_th * R_d<real_t>()) 
        / p_1000<real_t>() * (real_t(1) + r / eps<real_t>()),
      1 / (1 - R_over_c_p(r))
    )); 
  }

  // dtheta^star_drv from First Law for theta^star
  // TODO: document the derivation!
  phc_declare_funct_macro quantity<si::temperature, real_t> dtheta_drv(
    const quantity<si::temperature, real_t> T,
    const quantity<si::pressure, real_t> p,
    const quantity<mixing_ratio, real_t> r,
    const quantity<multiply_typeof_helper<si::mass_density, si::temperature>::type, real_t> rhod_th,
    const quantity<si::mass_density, real_t> rhod
  )
  {
    return - rhod_th / rhod * ( 
      // the 'liquid water' term
      l_v<real_t>(T) 
        / real_t(pow(1 + r, 2)) 
        / c_p(r) 
        / T 
/*
      +   
      // the 'virtual' term (TODO: as an option!)
      log(p / p_1000<real_t>()) 
        * R_d_over_c_pd<real_t>() 
        * (real_t(1) / eps<real_t>() - real_t(1) / ups<real_t>()) 
        * real_t(pow(real_t(1) + r / ups<real_t>(),-2)) 
*/
    );  
  }

};
