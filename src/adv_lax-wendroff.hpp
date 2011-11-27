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

template <class unit, typename real_t> 
class adv_lax_wendroff : public adv<unit, real_t> 
{
  public: const int stencil_extent() { return 3; }
  public: const int time_levels() { return 3; } // ?
  public: const int num_steps() { return 1; } // ?

  public: void op(Array<quantity<unit, real_t>, 3>* psi[], 
    const Range &i, 
    const Range &j, 
    const Range &k, 
    const int n, const int step,
    const Array<quantity<si::dimensionless, real_t>, 3> &Cx, 
    const Array<quantity<si::dimensionless, real_t>, 3> &, 
    const Array<quantity<si::dimensionless, real_t>, 3> &
  )
  {
    assert(step == 1);
    error_macro("not implemented yet")
/*
    (*psi[n+1])(i) -= 
      .5 * 
      (Cx(i+grid.p_half) + Cx(i-grid.m_half)) 
      * ( (*psi[n])(i+1) - (*psi[n])(i-1) );
*/
  }
};
#endif
