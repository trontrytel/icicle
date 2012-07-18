/** @file
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date July 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#pragma once
#include "phc.hpp"

namespace phc
{

  // molar mass of possible chemical compounds of the droplet
  // (to calculate from moles to mass)

  phc_declare_const_macro(M_SO2,  64.*1e-3, si::kilograms/si::moles) 
  phc_declare_const_macro(M_H2O2, 34.*1e-3, si::kilograms/si::moles) 
  phc_declare_const_macro(M_O3,   48.*1e-3, si::kilograms/si::moles)

  phc_declare_const_macro(M_HSO3, 81.*1e-3, si::kilograms/si::moles) 
  phc_declare_const_macro(M_SO3,  80.*1e-3, si::kilograms/si::moles) 
  phc_declare_const_macro(M_HSO4, 97.*1e-3, si::kilograms/si::moles)
  phc_declare_const_macro(M_SO4,  96.*1e-3, si::kilograms/si::moles)

};
