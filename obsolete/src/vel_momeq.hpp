/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011 - February 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#pragma once
#include "vel.hpp" 

template <typename real_t>
class vel_momeq : public vel<real_t>
{
  public: bool is_constant() const { return false; } 
};
