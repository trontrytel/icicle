/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#include "cfg.hpp"
#include "slv_parallel_threads.hpp"

#ifdef USE_BOOST_THREAD

template <typename real_t>
slv_parallel_threads<real_t>::slv_parallel_threads(
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
        nsd), nsd(nsd)
{
  int ncpu = boost::thread::hardware_concurrency();
  if (nsd > ncpu) warning_macro("using more threads (" << nsd << ") than CPUs/cores (" << ncpu << ")")
  b.reset(new boost::barrier(nsd));
}

template <typename real_t>
void slv_parallel_threads<real_t>::barrier()
{
  b->wait();
}

template <typename real_t>
void slv_parallel_threads<real_t>::integ_loop()
{
  boost::thread_group threads;
  for (int sd = 0; sd < nsd; ++sd) threads.add_thread(new boost::thread(
    &slv_parallel<real_t>::integ_loop_sd, this,sd
  ));
  threads.join_all();
}
  
#endif

// explicit instantiations
#if defined(USE_FLOAT)
template class slv_parallel_threads<float>;
#endif
#if defined(USE_DOUBLE)
template class slv_parallel_threads<double>;
#endif
#if defined(USE_LDOUBLE)
template class slv_parallel_threads<long double>;
#endif
