/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date January 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef TMP_HPP
#  define TMP_HPP

#  include "cmn.hpp" 
#  include "mtx.hpp"
#  include "grd.hpp"

template <typename real_t>
class tmp 
{
  private: ptr_vector<mtx::arr<real_t>> sclr_guard, vctr_guard;
  public: mtx::arr<real_t> **sclr, **vctr;
 
  public: tmp(
    int n_vctr, int n_sclr, const grd<real_t> &grid, int halomax,
    int i_min, int i_max, // TODO: encapsulate this info in some class?
    int j_min, int j_max, // TODO: ...
    int k_min, int k_max  // TODO: ...
  )
  {
    for (int n=0; n < n_vctr; ++n)
    {
      vctr_guard.push_back(new mtx::arr<real_t>(mtx::idx_ijk(
        // 3 x rng_vctr guarantees the temp space may be used for all dimensions
        grid.rng_vctr(i_min, i_max, halomax), 
        grid.rng_vctr(j_min, j_max, halomax), 
        grid.rng_vctr(k_min, k_max, halomax)  
      )));
    }
    vctr = vctr_guard.c_array();

    for (int n=0; n < n_sclr; ++n)
    {
      sclr_guard.push_back(new mtx::arr<real_t>(
        grid.rng_sclr(i_min, i_max, j_min, j_max, k_min, k_max, halomax)
      ));
    }
    sclr = sclr_guard.c_array();
  }
};
#endif
