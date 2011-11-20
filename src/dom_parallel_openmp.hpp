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
    out<unit, real_t> *output, vel<real_t> *velocity,
    int nx, quantity<si::length, real_t> dx, 
    int ny, quantity<si::length, real_t> dy, 
    int nz, quantity<si::length, real_t> dz, 
    quantity<si::time, real_t> dt,
    int nsd)
    : dom_parallel<unit, real_t>(fllbck, advsch, output, velocity, nx, dx, ny, dy, nz, dz, dt, nsd)
  {
    int ncpu = omp_get_num_procs();
    if (nsd > ncpu) warning_macro("using more threads (" << nsd << ") than CPUs/cores (" << ncpu << ")")
    omp_set_num_threads(nsd);
  }

  public: void barrier()
  {
#    pragma omp barrier
  }

  public: void integ_loop(unsigned long nt, quantity<si::time, real_t> dt)
  {
    int sd; 
#    pragma omp parallel private(sd)
    {   
      sd = omp_get_thread_num();
      this->integ_loop_sd(nt, dt, sd); 
    } 
  }
  
};
#  endif
#endif
