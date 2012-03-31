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

#  include "cmn.hpp"

// not using due to missing template<real_t> logic
//#  include <boost/units/systems/si/codata/physico-chemical_constants.hpp>

#  define declare_const_macro(name, value, unit) template <typename real_t> \
  static constexpr auto name() -> decltype(real_t(value) * unit) { return real_t(value) * unit; }

#  define derived_const_macro(name, value) template <typename real_t> \
  static constexpr auto name() -> decltype(value) { return value; }

namespace phc
{
  // acceleration due to gravity
  declare_const_macro(g, 9.81, si::metres_per_second_squared)

  // specific heat capacity of dry air
  declare_const_macro(c_pd, 1005, si::joules / si::kilograms / si::kelvins)

  // pressure in the definition of potential temperature
  declare_const_macro(p_0, 100000, si::pascals)

  // molar mass of dry air
  declare_const_macro(M_d, 0.02896, si::kilograms / si::moles)

  // universal gas constant
  declare_const_macro(R, 8.314472, si::joules / si::kelvins / si::moles)

  // gas constant for dry air
  derived_const_macro(R_d, R<real_t>() / M_d<real_t>())

  // Exner function exponent
  derived_const_macro(R_d_over_c_pd, R_d<real_t>() / c_pd<real_t>())
};

#  undef decl_const
#endif
