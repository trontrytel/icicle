/** @file
 *  @author Sylwester Arabas (University of Warsaw) <slayoo@igf.fuw.edu.pl>
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef ADV_LEAPFROG_HPP
#  define ADV_LEAPFROG_HPP

#  include "adv.hpp"

template <class unit, typename real_t> 
class adv_leapfrog : public adv<unit, real_t> 
{
  public: const int stencil_extent() { return 3; }
  public: const int time_levels() { return 3; }
  public: const int num_steps() { return 1; }
 
  public: void op_1D(Array<quantity<unit, real_t>, 3>* psi[], 
    const Range &i, const int n, const int step,
    const Array<quantity<si::dimensionless, real_t>, 3> &Cx, 
    const Array<quantity<si::dimensionless, real_t>, 3> &Cy, 
    const Array<quantity<si::dimensionless, real_t>, 3> &Cz,
    grd<real_t> &grid
  )
  {
    assert(step == 1);
    // TODO: is the averaging done here correct???
    (*psi[n+1])(i) -= 
      .5 * (Cx(i+grid.p_half) + Cx(i-grid.m_half)) // average Courant number
      * ( (*psi[n])(i+1) - (*psi[n])(i-1) );
  }
};
#endif
