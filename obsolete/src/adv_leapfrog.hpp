/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @brief definition of the adv_leapfrog class with an implementation of the leapfrog advection scheme
 */
#pragma once

#include "grd.hpp"
#include "adv_upstream.hpp"

#include <memory>
using std::unique_ptr;

/// @brief an implementation of the leapfrog advection scheme
template <typename real_t> 
class adv_leapfrog : public adv<real_t> 
{
  public: const int stencil_extent() const { return 3; }
  public: const int time_levels() const { return 3; }
 
  public: const real_t courant_max() const { return 1.; }
  public: const real_t courant_min() const { return .5; } ; //TODO Wicker Skamarock 2002 
 
  private: unique_ptr<adv_upstream<real_t>> fallback;
  public: adv_leapfrog();

  public: class op3D;
  typename adv<real_t>::op3D *factory(
    const mtx::idx &ijk,
    mtx::arr<real_t> **, 
    mtx::arr<real_t> **, 
    bool 
  ) const;
};
