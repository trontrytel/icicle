/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#pragma once
#include "slv_parallel_distmem.hpp"
#include "slv_parallel_serial.hpp"

#ifdef USE_BOOST_MPI
#  include <boost/mpi.hpp>
namespace mpi = boost::mpi;
#endif

template <typename real_t, class shrdmem_class>
class slv_parallel_distmem_mpi : public slv_parallel_distmem<real_t, shrdmem_class>
{
#ifdef USE_BOOST_MPI
  private: enum msgtype { msg_halo };
  private: mpi::environment *env; // perhaps MPI::Init_thread(MPI_THREAD_MULTIPLE) instead to support MPI+threads?
  private: mpi::communicator *comm;
  private: int mpi_init();
  public: ~slv_parallel_distmem_mpi();
  private: void sndrcv(int peer, int cnt, real_t *ibuf, real_t *obuf);
#endif

  public: slv_parallel_distmem_mpi(
    const stp<real_t> &setup, 
    out<real_t> &output, 
    int nsd
  ) 
#ifdef USE_BOOST_MPI
  ;
#else
  : slv_parallel_distmem<real_t, shrdmem_class>(setup, output, nsd, 0)
  {
    error_macro("please recompile icicle with -DMPI")
  }
#endif

  private: void distmem_barrier()
#ifdef USE_BOOST_MPI
  ;
#else
  {}
#endif
};
