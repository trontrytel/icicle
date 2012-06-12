/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#pragma once
#include "slv_parallel.hpp"
#ifdef USE_BOOST_THREAD
#  include <boost/thread.hpp>
#endif

template <typename real_t>
class slv_parallel_threads : public slv_parallel<real_t>
{
#ifdef USE_BOOST_THREAD
  private: unique_ptr<boost::barrier> b; 
#endif
  private: int nsd;

  public: slv_parallel_threads(
    const stp<real_t> &setup,
    out<real_t> &output,
    int i_min, int i_max,  
    int j_min, int j_max,  
    int k_min, int k_max,  
    int nsd)
#ifdef USE_BOOST_THREAD
  ;
#else
    : slv_parallel<real_t>(setup, output,
        i_min, i_max, 
        j_min, j_max,
        k_min, k_max, 
        nsd), nsd(nsd)
  {
    error_macro("please recompile icicle with -DTHREAD")
  }
#endif

  public: void barrier()
#ifdef USE_BOOST_THREAD
  ;
#else
  {}
#endif

  public: void integ_loop()
#ifdef USE_BOOST_THREAD
  ;
#else
  {}
#endif
};
