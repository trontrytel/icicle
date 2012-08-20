/** @file
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date August 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#pragma once
#include "phc.hpp"
#include "phc_maxwell-mason.hpp"

// parametrisations for 2-moment warm rain microphysics scheme
namespace phc
{
  namespace mm
  {
    // relaxation time for condensation/evaporation term for cloud and rain droplets
    // see Khvorostyaov at al 2001 eq. 5 
    phc_declare_funct_macro quantity<si::time, real_t> tau_relax(
      quantity<si::length, real_t> r, 
      quantity<divide_typeof_helper<si::dimensionless, si::volume>::type, real_t> N
    )
    {
      quantity<divide_typeof_helper<si::area, si::time>::type, real_t> D = phc::D<real_t>();
      return real_t(1) / (real_t(4 * M_PI) * D * N * r);
    }

    // d(es)/d(T) (Claus-Clapeyron eq. for ordinary atmospheric conditions + ideal gas assumption)
    // see Rogers & Yau: Claus-Clapeyron equation
    phc_declare_funct_macro quantity<divide_typeof_helper<si::pressure, si::temperature>::type, real_t> desdT(
      quantity<si::pressure, real_t> es,
      quantity<si::temperature, real_t> T 
    )
    {
    quantity<divide_typeof_helper<si::energy, si::mass>::type , real_t> lv = phc::l_v<real_t>(T);
    return lv * es / phc::R_v<real_t>() / T / T;
    }
  }
};
