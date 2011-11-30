/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef SLV_PARALLEL_SERIAL_HPP
#  define SLV_PARALLEL_SERIAL_HPP

#  include "slv_parallel.hpp"

template <class unit, typename real_t>
class slv_parallel_serial : public slv_parallel<unit, real_t>
{
  public: slv_parallel_serial(adv<unit, real_t> *fllbck, adv<unit, real_t> *advsch, 
    out<unit, real_t> *output, vel<real_t> *velocity, ini<real_t> *intcond,
    int i_min, int i_max, int nx, 
    int j_min, int j_max, int ny, 
    int k_min, int k_max, int nz, 
    grd<real_t> *grid, quantity<si::time, real_t> dt, int nsd
  )
    : slv_parallel<unit, real_t>(fllbck, advsch, output, velocity, intcond,
      i_min, i_max, nx, j_min, j_max, ny, k_min, k_max, nz, grid, dt, nsd
    )
  {
    assert(nsd == 1);
  }

  public: void barrier() 
  { 
  }

  public: void integ_loop(unsigned long nt, quantity<si::time, real_t> dt)
  {
    this->integ_loop_sd(nt, dt, 0); 
  }
};
#endif
