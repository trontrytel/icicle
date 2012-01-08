/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 *  @section DESCRIPTION
 *    an implementation of the Lax-Wendroff scheme following eq. 4 in
 *    Fischer and Renner 1972, MWR 100, no. 9
 *    (with the assumption that velocity field does not change with time)
 */
#ifndef ADV_LAX_WENDROFF_HPP
#  define ADV_LAX_WENDROFF_HPP

#  include "adv.hpp"

template <typename real_t> 
class adv_lax_wendroff : public adv<real_t> 
{
  adv_hack_macro // workaround for virtual template methods
  public: const int stencil_extent() { return 3; }
  public: const int time_levels() { return 3; } // ?
  public: const int num_steps() { return 1; } // ?

  public: 
  template <class idx>
  void op(int dim,
    arr<real_t>* psi[], 
    arr<real_t>* [], 
    arr<real_t>* [], 
    const Range &i, 
    const Range &j, 
    const Range &k, 
    const int n, const int step,
    const arr<real_t> * const Cx, 
    const arr<real_t> * const, 
    const arr<real_t> * const
  )
  {
    assert(step == 1);
    error_macro("not implemented yet")
  }
};
#endif
