/** @file
 *  @author Sylwester Arabas (University of Warsaw) <slayoo@igf.fuw.edu.pl>
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef DOM_PARALLEL_THREADS_HPP
#  define DOM_PARALLEL_THREADS_HPP

#  include "dom_parallel.hpp"
#  include <boost/thread.hpp>

template <class unit, typename real_t>
class dom_parallel_threads : public dom_parallel<unit, real_t>
{
  private: boost::barrier *b; 
  private: int nsd;

  public: dom_parallel_threads(adv<unit, real_t> *fllbck, adv<unit, real_t> *advsch, 
    out<unit, real_t> *output, int nx, int ny, int nz, int nsd)
    : dom_parallel<unit, real_t>(fllbck, advsch, output, nx, ny, nz, nsd), nsd(nsd)
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

  public: void integ_loop(unsigned long nt, 
    quantity<si::dimensionless, real_t> &Cx,
    quantity<si::dimensionless, real_t> &Cy,
    quantity<si::dimensionless, real_t> &Cz
  )
  {
    boost::thread_group threads;
    for (int sd = 0; sd < nsd; ++sd) threads.add_thread(new boost::thread(
      &dom_parallel<unit, real_t>::integ_loop_sd, this, nt, Cx, Cy, Cz, sd
    ));
    threads.join_all();
  }
  
};
#endif
