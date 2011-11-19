/** @file
 *  @author Sylwester Arabas (University of Warsaw) <slayoo@igf.fuw.edu.pl>
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef DOM_HPP
#  define DOM_HPP

#  include "config.hpp" // USE_* defines
#  include "common.hpp" // root class, error reporting
#  include "adv.hpp"

template <class unit, typename real_t>
class dom : root
{
  private: dom *left, *rght;
  public: void hook_neighbours(dom *l, dom *r) { left = l; rght = r; }

  private: virtual quantity<unit, real_t> data(int n, int i, int j, int k) = 0;

  public: virtual quantity<unit, real_t> left_nghbr_data(int n, int i, int j, int k) 
  { 
    return left->data(n, i, j, k); 
  }
  public: virtual quantity<unit, real_t> rght_nghbr_data(int n, int i, int j, int k) 
  { 
    return rght->data(n, i, j, k); 
  }

  public: bool choose_an(adv<unit, real_t> **a, int *n, int t, adv<unit, real_t> *advsch, adv<unit, real_t> *fllbck)
  {
    assert(advsch->time_levels() <= 3); // FIXME: support for other values
    bool fallback = (t == 0 && advsch->time_levels() == 3); 
    *a = fallback ? fllbck : advsch;
    *n = (*a)->time_levels() - 2;
    return fallback;
  }
  
  public: virtual void integ_loop(unsigned long nt, 
    quantity<si::dimensionless, real_t> &Cx,
    quantity<si::dimensionless, real_t> &Cy,
    quantity<si::dimensionless, real_t> &Cz
  ) = 0;
};

#endif
