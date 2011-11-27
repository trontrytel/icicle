/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef SLV_PARALLEL_THREADS_HPP
#  define SLV_PARALLEL_THREADS_HPP
#  ifdef USE_BOOST_THREAD

#    include "slv_parallel.hpp"
#    include <boost/thread.hpp>

template <class unit, typename real_t>
class slv_parallel_threads : public slv_parallel<unit, real_t>
{
  private: auto_ptr<boost::barrier> b; 
  private: int nsd;

  public: slv_parallel_threads(adv<unit, real_t> *fllbck, adv<unit, real_t> *advsch, 
    out<unit, real_t> *output, vel<real_t> *velocity, 
    int nx, int ny, int nz, 
    grd<real_t> *grid,
    quantity<si::time, real_t> dt,
    int nsd)
    : slv_parallel<unit, real_t>(fllbck, advsch, output, velocity, nx, ny, nz, grid, dt, nsd), nsd(nsd)
  {
    int ncpu = boost::thread::hardware_concurrency();
    if (nsd > ncpu) warning_macro("using more threads (" << nsd << ") than CPUs/cores (" << ncpu << ")")
    b.reset(new boost::barrier(nsd));
  }

  public: void barrier()
  {
    b->wait();
  }

  public: void integ_loop(unsigned long nt, quantity<si::time, real_t> dt)
  {
    boost::thread_group threads;
    for (int sd = 0; sd < nsd; ++sd) threads.add_thread(new boost::thread(
      &slv_parallel<unit, real_t>::integ_loop_sd, this, nt, dt, sd
    ));
    threads.join_all();
  }
  
};
#  endif
#endif
