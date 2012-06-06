/** @file
 *  @author Sylwester Arabas <sarabas@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date May 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef PHC_KAPPA_KOEHLER_HPP
#  define PHC_KAPPA_KOEHLER_HPP

#  include "phc.hpp"

namespace phc
{
  /// @brief equilibrium wet radius to the third power for a given:
  /// @arg dry radius to the third power
  /// @arg the solubility parameter kappa
  /// @arg ratio of abmient vapour density/pressure to saturation vapour density/pressure for pure water
  /// 
  /// the formula stems from applying the kappa-Koehler relation 
  /// (eq. 6 in @copydetails Petters_and_Kreidenweis_2007) to a stationarity
  /// condition for a vapour diffusion equation, which translates to
  /// zero differece of vapour density: rho_ambient - rho_surface = 0
  ///
  /// since rho_surface = rho_surface_pure_water * a(r_w, r_d, kappa)
  /// one can derive r_w as a function of r_d, kappa and the ratio
  /// of abmient and surface vapour densities
  ///
  /// for the kappa-Koehler parameterisation rw3 is linear with rd3
  phc_declare_funct_macro quantity<si::volume, real_t> rw3_eq(
    quantity<si::volume, real_t> rd3, 
    quantity<si::dimensionless, real_t> kappa,
    quantity<si::dimensionless, real_t> vap_ratio
  )
  {
    return rd3 * (1 - vap_ratio * (1 - kappa)) / (1 - vap_ratio);
  }
};

#endif
