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
class tmp : root
{
  private: unique_ptr<mtx::arr<real_t> > *sclr_guard;
  public: mtx::arr<real_t> **sclr;

  private: unique_ptr<mtx::arr<real_t> > *vctr_guard;
  public: mtx::arr<real_t> **vctr;
 
  public: tmp(int n_vctr, int n_sclr, grd<real_t> *grid, int halomax,
    int i_min, int i_max, // TODO: encapsulate this info in some class?
    int j_min, int j_max, // TODO: ...
    int k_min, int k_max  // TODO: ...
  )
  {
    vctr_guard = new unique_ptr<mtx::arr<real_t> >[n_vctr];
    vctr = new mtx::arr<real_t>*[n_vctr];
    for (int n=0; n < n_vctr; ++n)
    {
      vctr_guard[n].reset(new mtx::arr<real_t>(mtx::idx_ijk(
        // 3 x rng_vctr guarantees the temp space may be used for all dimensions
        grid->rng_vctr(i_min, i_max, halomax), 
        grid->rng_vctr(j_min, j_max, halomax), 
        grid->rng_vctr(k_min, k_max, halomax)  
      )));
      vctr[n] = vctr_guard[n].get();
    }

    sclr_guard = new unique_ptr<mtx::arr<real_t> >[n_sclr];
    sclr = new mtx::arr<real_t>*[n_sclr];
    for (int n=0; n < n_sclr; ++n)
    {
      sclr_guard[n].reset(new mtx::arr<real_t>(
        grid->rng_sclr(i_min, i_max, j_min, j_max, k_min, k_max, halomax)
      ));
      sclr[n] = sclr_guard[n].get();
    }
  }

  public: ~tmp() // TODO: get rid of it by using ptr_vector as in slv_serial.hpp
  {
    delete[] vctr;
    delete[] vctr_guard;
    delete[] sclr;
    delete[] sclr_guard;
  }
};
#endif
