/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef SLV_PARALLEL_OPENMP_HPP
#  define SLV_PARALLEL_OPENMP_HPP
#  ifdef _OPENMP

#    include "slv_parallel.hpp"
#    include <omp.h>

template <typename real_t>
class slv_parallel_openmp : public slv_parallel<real_t>
{
  public: slv_parallel_openmp(
    const stp<real_t> &setup, 
    out<real_t> &output,
    int i_min, int i_max, 
    int j_min, int j_max,  
    int k_min, int k_max,  
    int nsd)
    : slv_parallel<real_t>(setup, output,
      i_min, i_max, 
      j_min, j_max, 
      k_min, k_max, 
      nsd
  )
  {
    int ncpu = omp_get_num_procs();
    if (nsd > ncpu) warning_macro("using more threads (" << nsd << ") than CPUs/cores (" << ncpu << ")")
    omp_set_num_threads(nsd);
  }

  public: void barrier()
  {
#    pragma omp barrier
  }

  public: void integ_loop()
  {
    int sd; 
#    pragma omp parallel private(sd)
    {   
      sd = omp_get_thread_num();
      this->integ_loop_sd(sd); 
    } 
  }
  
};
#  endif
#endif
