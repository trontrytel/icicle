/** @file
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date August 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#pragma once
#include "phc.hpp"

namespace phc
{
  namespace kelvin
  {
    phc_declare_const_macro(sg_surf, 0.072, si::newton/si::metres);

    phc_declare_funct_macro quantity<si::length, real_t> A(
      quantity<si::temperature, real_t> T
    )
    {
      return real_t(2.) * sg_surf<real_t>() / phc::R_v<real_t>() / T / phc::rho_w<real_t>();
    }
  }
};
