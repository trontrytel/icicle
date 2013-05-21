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
class slv_parallel_serial : public slv_parallel<real_t>
{
  public: slv_parallel_serial(
    const stp<real_t> &setup, 
    out<real_t> &output,
    int i_min, int i_max, 
    int j_min, int j_max, 
    int k_min, int k_max, 
    int nsd = 1
  );

  public: void barrier() {}

  public: void integ_loop();
};
