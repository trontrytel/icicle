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
  namespace react
  {
  // reactivity constants for S(IV) -> S(VI) chemical reactions
  // Seinfeld & Pandis

    phc_declare_const_macro(R_S_H2O2_k, 7.5  * 1e7  *1e-6, si::cubic_metres*si::cubic_metres/si::moles/si::moles/si::seconds) 
    phc_declare_const_macro(R_S_H2O2_K, 13.         *1e-6, si::cubic_metres/si::moles) 

    phc_declare_const_macro(R_S_O3_k0,   2.4 * 1e4  *1e-3, si::cubic_metres/si::moles/si::seconds) 
    phc_declare_const_macro(R_S_O3_k1,   3.7 * 1e5  *1e-3, si::cubic_metres/si::moles/si::seconds) 
    phc_declare_const_macro(R_S_O3_k2,   1.5 * 1e9  *1e-3, si::cubic_metres/si::moles/si::seconds) 
  }
};
