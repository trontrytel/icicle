/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef SLV_HPP
#  define SLV_HPP

#  include "cmn.hpp" // root class, error reporting
#  include "adv.hpp"

template <typename real_t>
class slv : root
{
  // left->i_min, rght->i_max, fore->j_min, hind->j_max, base->k_min, apex->k_max (enums may not be ++)
  public: static const int first=0, left=first, rght=1, fore=2, hind=3, base=4, apex=5, last=apex;
  private: slv *nghbrs[6];

  public: virtual void hook_neighbour(int s, slv *n) 
  { 
    nghbrs[s] = n; 
  }

  public: virtual typename arr<real_t>::arr_ret data(int n, const idx &idx) = 0;

  public: typename arr<real_t>::arr_ret nghbr_data(int side, int n, const idx &idx)
  {
    nghbrs[side]->sync(n);
    return nghbrs[side]->data(n, idx);
  }

  public: bool choose_an(adv<real_t> **a, int *n, int t, 
    adv<real_t> *advsch, adv<real_t> *fllbck
  )
  {
    assert(advsch->time_levels() <= 3); // FIXME: support for other values
    bool fallback = (t == 0 && advsch->time_levels() == 3); 
    *a = fallback ? fllbck : advsch;
    *n = (*a)->time_levels() - 2;
    return fallback;
  }

  public: virtual void integ_loop() = 0;
  public: virtual void sync(int n) {};
};

#endif
