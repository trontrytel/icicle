/** @file
 *  @author Sylwester Arabas (University of Warsaw) <slayoo@igf.fuw.edu.pl>
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef DOM_PARALLEL_THREADS_HPP
#  define DOM_PARALLEL_THREADS_HPP
#  ifdef USE_BOOST_THREAD

#    include "dom_parallel.hpp"
#    include <boost/thread.hpp>

template <class unit, typename real_t>
class dom_parallel_threads : public dom_parallel<unit, real_t>
{
  private: boost::barrier *b; 
  private: int nsd;

  public: dom_parallel_threads(adv<unit, real_t> *fllbck, adv<unit, real_t> *advsch, 
    out<unit, real_t> *output, vel<real_t> *velocity, 
    int nx, int ny, int nz, 
    grd<real_t> *grid,
    quantity<si::time, real_t> dt,
    int nsd)
    : dom_parallel<unit, real_t>(fllbck, advsch, output, velocity, nx, ny, nz, grid, dt, nsd), nsd(nsd)
  {
    int ncpu = boost::thread::hardware_concurrency();
    if (nsd > ncpu) warning_macro("using more threads (" << nsd << ") than CPUs/cores (" << ncpu << ")")
    b = new boost::barrier(nsd);
  }

  public: ~dom_parallel_threads()
  {
    delete b;
  }

  public: void barrier()
  {
    b->wait();
  }

  public: void integ_loop(unsigned long nt, quantity<si::time, real_t> dt)
  {
    boost::thread_group threads;
    for (int sd = 0; sd < nsd; ++sd) threads.add_thread(new boost::thread(
      &dom_parallel<unit, real_t>::integ_loop_sd, this, nt, dt, sd
    ));
    threads.join_all();
  }
  
};
#  endif
#endif
