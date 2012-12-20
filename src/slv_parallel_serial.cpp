/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#include "slv_parallel_serial.hpp"

template <typename real_t>
slv_parallel_serial<real_t>::slv_parallel_serial(
  const stp<real_t> &setup, 
  out<real_t> &output,
  int i_min, int i_max, 
  int j_min, int j_max, 
  int k_min, int k_max, 
  int nsd
)
: slv_parallel<real_t>(setup, output, i_min, i_max, j_min, j_max, k_min, k_max, 1)
{ 
  assert(nsd == 1);
}

template <typename real_t>
void slv_parallel_serial<real_t>::integ_loop()
{
  this->integ_loop_sd(0); 
}

// explicit instantiations
#define ICICLE_INSTANTIATE_CLASS slv_parallel_serial
#include "cmn/cmn_instant.hpp"
