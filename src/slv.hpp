/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef SLV_HPP
#  define SLV_HPP

#  include "config.hpp" // USE_* defines
#  include "common.hpp" // root class, error reporting
#  include "adv.hpp"

template <typename real_t>
class slv : root
{
  // left->i_min, rght->i_max, fore->j_min, hind->j_max, base->k_min, apex->k_max
  public: enum side { first, left=first, rght, fore, hind, base, apex, last=apex };
  private: slv *nghbrs[6];

  public: virtual void hook_neighbour(side s, slv *n) 
  { 
    nghbrs[s] = n; 
  }

  public: virtual Array<real_t, 3> data(
    int n, const Range &i, const Range &j, const Range &k
  ) = 0;

  public: Array<real_t, 3> nghbr_data(
    side s, int n, const Range &i, const Range &j, const Range &k
  )
  {
    nghbrs[s]->sync(n);
    return nghbrs[s]->data(n, i, j, k);
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

  public: virtual void integ_loop(unsigned long nt, quantity<si::time, real_t> dt) = 0;
  public: virtual void sync(int n) {};
};

#endif
