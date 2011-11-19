/** @file
 *  @author Sylwester Arabas (University of Warsaw) <slayoo@igf.fuw.edu.pl>
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef DOM_PARALLEL_DISTMEM_HPP
#  define DOM_PARALLEL_DISTMEM_HPP

#  include "dom_parallel.hpp"

template <class unit, typename real_t, class shrdmem_class>
class dom_parallel_openmp : public shrdmem_class
{
  public: dom_parallel_distmem(adv<unit, real_t> *fllbck, adv<unit, real_t> *advsch, 
    out<unit, real_t> *output, int nx, int ny, int nz, int nsd)
    : shrdmem_class(fllbck, advsch, output, nx, ny, nz, nsd)
  {
  }

  public: void barrier()
  {
  }

  public: void integ_loop(unsigned long nt, 
    quantity<si::dimensionless, real_t> &Cx,
    quantity<si::dimensionless, real_t> &Cy,
    quantity<si::dimensionless, real_t> &Cz
  )
  {
  }
  
};
#endif
