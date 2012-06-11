/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#include "cfg.hpp"
#include "slv_parallel_openmp.hpp"
#ifdef _OPENMP
#  include <omp.h>

template <typename real_t>
slv_parallel_openmp<real_t>::slv_parallel_openmp(
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

template <typename real_t>
void slv_parallel_openmp<real_t>::barrier()
{
#  pragma omp barrier
}

template <typename real_t>
void slv_parallel_openmp<real_t>::integ_loop()
{
  int sd; 
#  pragma omp parallel private(sd)
  {   
    sd = omp_get_thread_num();
    this->integ_loop_sd(sd); 
  } 
}
#endif

// explicit instantiations
#if defined(USE_FLOAT)
template class slv_parallel_openmp<float>;
#endif
#if defined(USE_DOUBLE)
template class slv_parallel_openmp<double>;
#endif
#if defined(USE_LDOUBLE)
template class slv_parallel_openmp<long double>;
#endif
