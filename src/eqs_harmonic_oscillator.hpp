/** @file
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date March 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 @  @brief contains definition of eqs_harmonic_oscillator
 */
#pragma once

//#include "cmn.hpp"
#include "eqs.hpp"

/// @brief a minimalistic model of two coupled harmonic oscillators
///   (consult eq. 28 in Smolarkiewicz 2006, IJNMF)
template <typename real_t>
class eqs_harmonic_oscillator : public eqs<real_t> 
{
  private: class restoring_force;
  public: eqs_harmonic_oscillator(quantity<si::frequency, real_t> omega);
};
