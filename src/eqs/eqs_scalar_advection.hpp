/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date January 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @brief contains definition of the eqs_scalar_advection class - the simplest equation system
 */
#pragma once

#include "eqs.hpp"

/// @brief the simplest equation system consisting of a single homogeneous transport equation
template <typename real_t>
class eqs_scalar_advection : public eqs<real_t>
{
  public: eqs_scalar_advection();
};