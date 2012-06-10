/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011 - March 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @brief definition of the adv_mpdata_fct class with an implementation of the flux-corrected MPDATA scheme
 */
#pragma once

#include "adv_mpdata.hpp"

/** @brief
 *  Flux Corrected Transport (aka non-oscillatory, monotonic, sign-preserving)
 *  option for MPDATA (for the transport of positive scalars only)
 */
template <typename real_t> 
class adv_mpdata_fct : public adv_mpdata<real_t> 
{
  public: const int stencil_extent() const
  {
    int halo = (adv_mpdata<real_t>::stencil_extent() - 1) / 2; 
    halo += 1; // cf. ii, jj & kk below
    return 1 + 2 * halo;
  }
  public: const int num_vctr_caches() const
  { 
    return 1 + (iord < 3 ? 3 : 6); 
  }
  public: const int num_sclr_caches() const { return 2; } 

  private: int iord;
  private: bool cross_terms, third_order;
  public: adv_mpdata_fct(int iord, bool cross_terms, bool third_order);

  public: template <class mpdata_op3D_class> class op3D;
  
  public: typename adv<real_t>::op3D *factory(
    const mtx::idx &ijk,
    mtx::arr<real_t> **tmp_s,
    mtx::arr<real_t> **tmp_v,
    bool positive_definite
  ) const;
};
