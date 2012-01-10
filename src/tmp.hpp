/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date January 2012
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef TMP_HPP
#  define TMP_HPP

#  include "cmn.hpp" // root class, error reporting
#  include "arr.hpp"

template <typename real_t>
class tmp : root
{
  private: auto_ptr<arr<real_t> > *sclr_guard;
  public: arr<real_t> **sclr;

  private: auto_ptr<arr<real_t> > *vctr_guard;
  public: arr<real_t> **vctr;
 
  public: tmp(int n_vctr, int n_sclr, grd<real_t> *grid, int halo,
    int i_min, int i_max, // TODO: encapsulate this info in some class?
    int j_min, int j_max, // TODO: ...
    int k_min, int k_max  // TODO: ...
  )
  {
    vctr_guard = new auto_ptr<arr<real_t> >[n_vctr];
    vctr = new arr<real_t>*[n_vctr];
    for (int n=0; n < n_vctr; ++n)
    {
      vctr_guard[n].reset(new arr<real_t>(
        // 3 x rng_vctr guarantees the temp space may be used for all dimensions
        grid->rng_vctr(i_min, i_max, halo), 
        grid->rng_vctr(j_min, j_max, halo), 
        grid->rng_vctr(k_min, k_max, halo)  
      ));
      vctr[n] = vctr_guard[n].get();
    }

    sclr_guard = new auto_ptr<arr<real_t> >[n_sclr];
    sclr = new arr<real_t>*[n_sclr];
    for (int n=0; n < n_sclr; ++n)
    {
      sclr_guard[n].reset(new arr<real_t>(
        grid->rng_sclr(i_min, i_max, halo),
        grid->rng_sclr(j_min, j_max, halo),
        grid->rng_sclr(k_min, k_max, halo)
      ));
      sclr[n] = sclr_guard[n].get();
    }
  }

  public: ~tmp()
  {
    delete[] vctr;
    delete[] vctr_guard;
    delete[] sclr;
    delete[] sclr_guard;
  }
};
#endif