/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
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

  public: virtual typename mtx::arr<real_t>::type data(int e, int n, const mtx::idx &idx) = 0;

  public: typename mtx::arr<real_t>::type nghbr_data(int side, int e, int n, const mtx::idx &idx)
  {
    nghbrs[side]->sync(e, n);
    return nghbrs[side]->data(e, n, idx);
  }

  public: virtual void integ_loop() = 0;
  public: virtual void sync(int e, int n) {};
};

#endif
