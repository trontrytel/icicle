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
  namespace dissociation
  {
  // dissociation constants for chemical compounds of droplets

  // TODO dissociation constants depend on temperature 
  // (provided values are for T=298K)

    phc_declare_const_macro(K_H2O, 1e-16    *1e3, si::moles / si::cubic_metres) 
    phc_declare_const_macro(K_SO2, 1.3*1e-2 *1e3, si::moles / si::cubic_metres) 
    phc_declare_const_macro(K_HSO3,6.6*1e-8 *1e3, si::moles / si::cubic_metres)
  }
};
