/** @file
 *  @author Sylwester Arabas (University of Warsaw) <slayoo@igf.fuw.edu.pl>
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef DOM_PARALLEL_OPENMP_HPP
#  define DOM_PARALLEL_OPENMP_HPP
#  ifdef _OPENMP

#    include "dom_parallel.hpp"
#    include <omp.h>

template <class unit, typename real_t>
class dom_parallel_openmp : public dom_parallel<unit, real_t>
{
  public: dom_parallel_openmp(adv<unit, real_t> *fllbck, adv<unit, real_t> *advsch, 
    out<unit, real_t> *output, int nx, int ny, int nz, int nsd)
    : dom_parallel<unit, real_t>(fllbck, advsch, output, nx, ny, nz, nsd)
  {
    int ncpu = omp_get_num_procs();
    if (nsd > ncpu) warning_macro("using more threads (" << nsd << ") than CPUs/cores (" << ncpu << ")")
    omp_set_num_threads(nsd);
  }

  public: void barrier()
  {
#    pragma omp barrier
  }

  public: void integ_loop(unsigned long nt, 
    quantity<si::dimensionless, real_t> &Cx,
    quantity<si::dimensionless, real_t> &Cy,
    quantity<si::dimensionless, real_t> &Cz
  )
  {
    int sd; 
#    pragma omp parallel private(sd)
    {   
      sd = omp_get_thread_num();
      this->integ_loop_sd(nt, Cx, Cy, Cz, sd); 
    } 
  }
  
};
#  endif
#endif
