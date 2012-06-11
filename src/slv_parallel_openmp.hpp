/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#pragma once
#include "slv_parallel.hpp"

template <typename real_t>
class slv_parallel_openmp : public slv_parallel<real_t>
{
  public: slv_parallel_openmp(
    const stp<real_t> &setup, 
    out<real_t> &output,
    int i_min, int i_max, 
    int j_min, int j_max,  
    int k_min, int k_max,  
    int nsd
  )
#ifdef _OPENMP
  ;
#else
  : slv_parallel<real_t>(setup, output, i_min, i_max, j_min, j_max, k_min, k_max, nsd)
  {
    error_macro("please recompile icicle with -DOPENMP")
  }
#endif

  public: void barrier();

  public: void integ_loop();
};
