/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date March 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @brief contains definition of the phc class encapsulating 
 *    a catalogue of physical constants
 */
#ifndef PHC_HPP
#  define PHC_HPP

#  define decltype_return(expr) -> decltype(expr) { return expr; }

#  define declare_const_macro(name, value, unit) template <typename real_t> \
  static constexpr auto name() decltype_return(real_t(value) * unit)

#  define derived_const_macro(name, value) template <typename real_t> \
  static constexpr auto name() decltype_return(value)

#  define declare_funct_macro template <typename real_t> constexpr auto

// TODO: would changing namespace to class and hence marking all members defined below as inline help?
//       (but that excludes definition of namespace embers in multiple files)
namespace phc
{
  typedef si::dimensionless mixing_ratio;

  // acceleration due to gravity
  declare_const_macro(g, 9.81, si::metres_per_second_squared)

  // specific heat capacities
  declare_const_macro(c_pd, 1005, si::joules / si::kilograms / si::kelvins) // dry air
  declare_const_macro(c_pv, 1850, si::joules / si::kilograms / si::kelvins) // water vapour
  declare_const_macro(c_pw, 4218, si::joules / si::kilograms / si::kelvins) // liquid water

  // pressure in the definition of potential temperature
  declare_const_macro(p_1000, 100000, si::pascals)

  // molar masses
  declare_const_macro(M_d, 0.02896, si::kilograms / si::moles) // dry air
  declare_const_macro(M_v, 0.01802, si::kilograms / si::moles) // water vapour
  derived_const_macro(eps, M_v<real_t>() / M_d<real_t>())

  // universal gas constant (i.e. the Boltzmann times the Avogadro constants)
  declare_const_macro(kaBoNA, 8.314472, si::joules / si::kelvins / si::moles)

  // gas constants
  derived_const_macro(R_d, kaBoNA<real_t>() / M_d<real_t>()) // dry air
  derived_const_macro(R_v, kaBoNA<real_t>() / M_v<real_t>()) // water vapour

  // Exner function exponent
  derived_const_macro(R_d_over_c_pd, R_d<real_t>() / c_pd<real_t>())

  // water triple point parameters
  declare_const_macro(p_tri, 611.73, si::pascals) // pressure
  declare_const_macro(T_tri, 273.16, si::kelvins) // temperature
  declare_const_macro(l_tri, 2.5e6, si::joules / si::kilograms) // latent heat of evaporation

  // mixing rule for extensive quantitites (i.e. using mass mixing ratio)
  template <typename real_t, typename quant>
  auto constexpr mix(quant dry, quant vap, quantity<mixing_ratio, real_t> r)
    decltype_return((dry + r * vap) / (1 + r))

  // gas constant for moist air
  declare_funct_macro R(quantity<mixing_ratio, real_t> r)
    decltype_return(mix(R_d<real_t>(), R_v<real_t>(), r))
 
  // specific heat capacity of moist air
  declare_funct_macro c_p(quantity<mixing_ratio, real_t> r)
    decltype_return(mix(c_pd<real_t>(), c_pv<real_t>(), r))

  // saturation vapour pressure for water assuming constant c_p_v and c_p_w 
  // with constants taken at triple point
  // (solution to the Clausius-Clapeyron equation assuming rho_vapour << rho_liquid)
  declare_funct_macro p_vs(quantity<si::temperature, real_t> T)
    decltype_return(p_tri<real_t>() * exp(
      (l_tri<real_t>() + (c_pw<real_t>() - c_pv<real_t>()) * T_tri<real_t>()) / R_v<real_t>() * (1 / T_tri<real_t>() - 1 / T)
      - (c_pw<real_t>() - c_pv<real_t>()) / R_v<real_t>() * log(T / T_tri<real_t>())
    ))

  // saturation vapour mixing ratio for water
  declare_funct_macro r_vs(quantity<si::temperature, real_t> T, quantity<si::pressure, real_t> p)
    decltype_return(eps<real_t>() / (p / p_vs<real_t>(T) - 1))

  // Exner function for dry air
  declare_funct_macro exner(quantity<si::pressure, real_t> p)
    decltype_return(pow(p / p_1000<real_t>(), R_d_over_c_pd<real_t>()))

  // Exner function for moist air
  declare_funct_macro exner(quantity<si::pressure, real_t> p, quantity<mixing_ratio, real_t> r)
    decltype_return(pow(p / p_1000<real_t>(), R<real_t>(r) / c_p<real_t>(r)))

  // water vapour partial pressure as a function of mixing ratio
  declare_funct_macro p_v(quantity<si::pressure, real_t> p, quantity<mixing_ratio, real_t> r)
    decltype_return(p * r / (r + eps<real_t>()))
};

#  undef decltype_return
#  undef declare_funct_macro
#  undef declare_const_macro
#  undef derived_const_macro
#endif
