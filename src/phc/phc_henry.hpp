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
  namespace henry 
  {
  // Henry law coefficients for dissolving gaseous phase into water droplets

  //TODO Henrys coeffs should depend on temperature
  //TODO Henrys coeffs should depend on dissociation
  //TODO ude boost units for calculating to si

  //traditionally(?) Henrys constants are given in moles/liter/standard_atmosphere 
  //hence the need of /1.01325*1e3*1e-5 
    phc_declare_const_macro(H_SO2, 1.23      /1.01325*1e3*1e-5, si::mole/si::cubic_metres/si::pascals) 
    phc_declare_const_macro(H_H2O2,7.45*1e4  /1.01325*1e3*1e-5, si::moles/si::cubic_metres/si::pascals) 
    phc_declare_const_macro(H_O3,  1.13*1e-2 /1.01325*1e3*1e-5, si::moles/si::cubic_metres/si::pascals)
  }
};
