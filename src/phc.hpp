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

#  define phc_decltype_return(expr) -> decltype(expr) { return expr; }

#  define phc_declare_const_macro(name, value, unit) template <typename real_t> \
  static constexpr auto name() phc_decltype_return(real_t(value) * unit)

#  define phc_derived_const_macro(name, value) template <typename real_t> \
  static constexpr auto name() phc_decltype_return(value)

#  define phc_declare_funct_macro template <typename real_t> constexpr 

// TODO: the same functions with Blitz arguments - how to automate???

// TODO: would changing namespace to class and hence marking all members defined below as inline help?
//       (but that excludes definition of namespace embers in multiple files)
namespace phc
{
  typedef si::dimensionless mixing_ratio;

  // acceleration due to gravity
  phc_declare_const_macro(g, 9.81, si::metres_per_second_squared)

  // specific heat capacities
  phc_declare_const_macro(c_pd, 1005, si::joules / si::kilograms / si::kelvins) // dry air
  phc_declare_const_macro(c_pv, 1850, si::joules / si::kilograms / si::kelvins) // water vapour
  phc_declare_const_macro(c_pw, 4218, si::joules / si::kilograms / si::kelvins) // liquid water

  // pressure in the definition of potential temperature
  phc_declare_const_macro(p_1000, 100000, si::pascals)

  // molar masses
  phc_declare_const_macro(M_d, 0.02896, si::kilograms / si::moles) // dry air
  phc_declare_const_macro(M_v, 0.01802, si::kilograms / si::moles) // water vapour
  phc_derived_const_macro(eps, M_v<real_t>() / M_d<real_t>()) // aka epsilon
  phc_derived_const_macro(ups, c_pd<real_t>() / c_pv<real_t>()) // aka upsilon

  // universal gas constant (i.e. the Boltzmann times the Avogadro constants)
  phc_declare_const_macro(kaBoNA, 8.314472, si::joules / si::kelvins / si::moles)

  // gas constants
  phc_derived_const_macro(R_d, kaBoNA<real_t>() / M_d<real_t>()) // dry air
  phc_derived_const_macro(R_v, kaBoNA<real_t>() / M_v<real_t>()) // water vapour

  // Exner function exponent
  phc_derived_const_macro(R_d_over_c_pd, R_d<real_t>() / c_pd<real_t>())

  // water triple point parameters
  phc_declare_const_macro(p_tri, 611.73, si::pascals) // pressure
  phc_declare_const_macro(T_tri, 273.16, si::kelvins) // temperature
  phc_declare_const_macro(l_tri, 2.5e6, si::joules / si::kilograms) // latent heat of evaporation

  // mixing rule for extensive quantitites (i.e. using mass mixing ratio)
  template <typename real_t, typename quant>
  auto constexpr mix(quant dry, quant vap, quantity<mixing_ratio, real_t> r)
    phc_decltype_return((dry + r * vap) / (1 + r))

  // gas constant for moist air
  phc_declare_funct_macro auto R(quantity<mixing_ratio, real_t> r)
    phc_decltype_return(mix(R_d<real_t>(), R_v<real_t>(), r))
 
  // specific heat capacity of moist air
  phc_declare_funct_macro auto c_p(quantity<mixing_ratio, real_t> r)
    phc_decltype_return(mix(c_pd<real_t>(), c_pv<real_t>(), r))

  // latent heat for constant c_p
  phc_declare_funct_macro quantity<divide_typeof_helper<si::energy, si::mass>::type , real_t> l_v(quantity<si::temperature, real_t> T)
  {
    return l_tri<real_t>() + (c_pv<real_t>() - c_pw<real_t>()) * (T - T_tri<real_t>());
  }

  // saturation vapour pressure for water assuming constant c_p_v and c_p_w 
  // with constants taken at triple point
  // (solution to the Clausius-Clapeyron equation assuming rho_vapour << rho_liquid)
  phc_declare_funct_macro quantity<si::pressure, real_t> p_vs(quantity<si::temperature, real_t> T)
  {
    return p_tri<real_t>() * exp(
      (l_tri<real_t>() + (c_pw<real_t>() - c_pv<real_t>()) * T_tri<real_t>()) / R_v<real_t>() * (real_t(1) / T_tri<real_t>() - real_t(1) / T)
      - (c_pw<real_t>() - c_pv<real_t>()) / R_v<real_t>() * log(T / T_tri<real_t>())
    );
  }

  // saturation vapour mixing ratio for water
  phc_declare_funct_macro quantity<mixing_ratio, real_t> r_vs(quantity<si::temperature, real_t> T, quantity<si::pressure, real_t> p)
  {
    return eps<real_t>() / (p / p_vs<real_t>(T) - 1);
  }

  // Exner function for dry air
  phc_declare_funct_macro quantity<si::dimensionless, real_t> exner(quantity<si::pressure, real_t> p)
  {
    return pow(p / p_1000<real_t>(), R_d_over_c_pd<real_t>());
  }

  // Exner function exponent for moist air
  phc_declare_funct_macro quantity<si::dimensionless, real_t> R_over_c_p(quantity<mixing_ratio, real_t> r)
  {
    return R<real_t>(r) / c_p<real_t>(r);
  }

  // Exner function for moist air
  phc_declare_funct_macro quantity<si::dimensionless, real_t> exner(quantity<si::pressure, real_t> p, quantity<mixing_ratio, real_t> r)
  {
    return pow(p / p_1000<real_t>(), R_over_c_p<real_t>(r));
  }

  // water vapour partial pressure as a function of mixing ratio
  phc_declare_funct_macro quantity<si::pressure, real_t> p_v(quantity<si::pressure, real_t> p, quantity<mixing_ratio, real_t> r)
  {
    return p * r / (r + eps<real_t>());
  }
};

#endif
