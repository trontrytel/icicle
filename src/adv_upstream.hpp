/** @file
 *  @author Sylwester Arabas (University of Warsaw) <slayoo@igf.fuw.edu.pl>
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef ADV_UPSTREAM_HPP
#  define ADV_UPSTREAM_HPP

#  include "adv.hpp"

template <class unit, typename real_t>
class adv_upstream : public adv<unit, real_t>
{
  public: const int stencil_extent() { return 3; }
  public: const int time_levels() { return 2; }
  public: const int num_steps() { return 1; }

  public: void op_1D(Array<quantity<unit, real_t>, 3> *psi[], const Range &i,
    const int n, const quantity<si::dimensionless, real_t> &C, const int step)
  {
    assert(step == 1);
    // TODO: the same using where()!
    (*psi[n+1])(i) -= C * ( 
      C > 0.  
        ? (*psi[n])(i  ) - (*psi[n])(i-1) 
        : (*psi[n])(i+1) - (*psi[n])(i  )
    );  
  }
};
#endif
