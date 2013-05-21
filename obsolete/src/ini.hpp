/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @brief contains definition of the @ref ini class - a base class for initial conditions
 */
#pragma once
#include "grd.hpp"
#include "mtx.hpp"

#include <string>
using std::string;

/// @brief base class for initial conditions
template <typename real_t>
class ini 
{
  public: virtual void populate_scalar_field(
    const string &varname,
    const mtx::idx &ijk,
    mtx::arr<real_t> &data
  ) const = 0;
};
