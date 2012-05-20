/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date March 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @brief contains definitions of thermodynamic relations that result
 *    from the assumption of constant specific heat capacitise
 *    (i.e. formulas consistent with perfect gas model, but not 
 *    neccesarily the most efficient ones)
 */
#ifndef PHC_CONST_CP_HPP
#  define PHC_CONST_CP_HPP

#  include "phc.hpp"

// TODO: perhaps a template parameter (svp?) that would allow
//       swithing between different formulations 
//       (e.g. const_sp, flatau_1992, wmo, etc...)

namespace phc
{
  // saturation vapour pressure for water assuming constant c_p_v and c_p_w 
  // with constants taken at triple point
  // (solution to the Clausius-Clapeyron equation assuming rho_vapour << rho_liquid)
  phc_declare_funct_macro quantity<si::pressure, real_t> p_vs(
    quantity<si::temperature, real_t> T
  )
  {
    return p_tri<real_t>() * exp(
      (l_tri<real_t>() + (c_pw<real_t>() - c_pv<real_t>()) * T_tri<real_t>()) / R_v<real_t>() * (real_t(1) / T_tri<real_t>() - real_t(1) / T)
      - (c_pw<real_t>() - c_pv<real_t>()) / R_v<real_t>() * log(T / T_tri<real_t>())
    );
  }

  // saturation vapour mixing ratio for water as a function of pressure and temperature
  phc_declare_funct_macro quantity<mixing_ratio, real_t> r_vs(
    quantity<si::temperature, real_t> T, 
    quantity<si::pressure, real_t> p
  )
  {
    return eps<real_t>() / (p / p_vs<real_t>(T) - 1);
  }

  // latent heat for constant c_p
  phc_declare_funct_macro quantity<divide_typeof_helper<si::energy, si::mass>::type , real_t> l_v(
    quantity<si::temperature, real_t> T
  )
  {
    return l_tri<real_t>() + (c_pv<real_t>() - c_pw<real_t>()) * (T - T_tri<real_t>());
  }
};

#endif
