/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#pragma once
//TODO: split into .hpp and .cpp
#  ifdef USE_BOOST_MPI

#  include "slv_parallel_distmem.hpp"

#  include <boost/mpi.hpp>
namespace mpi = boost::mpi;

template <typename real_t, class shrdmem_class>
class slv_parallel_distmem_mpi : public slv_parallel_distmem<real_t, shrdmem_class>
{
  private: enum msgtype { msg_halo };
  private: mpi::environment *env; // perhaps MPI::Init_thread(MPI_THREAD_MULTIPLE) instead to support MPI+threads?
  private: mpi::communicator *comm;

  /// MPI init for using rank as an argument to the parent class c-tor
  private: int mpi_init()
  {
    env = new mpi::environment(); 
    comm = new mpi::communicator();
    return comm->rank();
  }

  public: slv_parallel_distmem_mpi(
    const stp<real_t> &setup, 
    out<real_t> &output, 
    int nsd
  ) : slv_parallel_distmem<real_t, shrdmem_class>(setup, output, nsd, mpi_init())
  {
    if (nsd != comm->size()) 
      error_macro("MPI is using different number of nodes (" << comm->size() << ") than specified (" << nsd << ")")
  }

  public: ~slv_parallel_distmem_mpi()
  {
    delete comm;
    delete env; 
  }

  private: void distmem_barrier()
  {
    comm->barrier();
  }

  private: void sndrcv(int peer, int cnt, real_t *ibuf, real_t *obuf)
  {
    mpi::request reqs[2];
    int r = 0;
    // isend() / irecv() are non-blocking
    reqs[r++] = comm->isend(peer, msg_halo, obuf, cnt); 
    reqs[r++] = comm->irecv(peer, msg_halo, ibuf, cnt);
    mpi::wait_all(reqs, reqs + r);
  }
};
#  endif
