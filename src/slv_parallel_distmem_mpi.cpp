/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#include "slv_parallel_distmem_mpi.hpp"
#if defined(USE_BOOST_MPI)

/// MPI init for using rank as an argument to the parent class c-tor
template <typename real_t, class shrdmem_class>
int slv_parallel_distmem_mpi<real_t, shrdmem_class>::mpi_init()
{
  env = new mpi::environment(); 
  comm = new mpi::communicator();
  return comm->rank();
}

template <typename real_t, class shrdmem_class>
slv_parallel_distmem_mpi<real_t, shrdmem_class>::slv_parallel_distmem_mpi(
  const stp<real_t> &setup, 
  out<real_t> &output, 
  int nsd
) : slv_parallel_distmem<real_t, shrdmem_class>(setup, output, nsd, mpi_init())
{
  if (nsd != comm->size()) 
    error_macro("MPI is using different number of nodes (" << comm->size() << ") than specified (" << nsd << ")")
}

template <typename real_t, class shrdmem_class>
slv_parallel_distmem_mpi<real_t, shrdmem_class>::~slv_parallel_distmem_mpi()
{
  delete comm;
  delete env; 
}

template <typename real_t, class shrdmem_class>
void slv_parallel_distmem_mpi<real_t, shrdmem_class>::distmem_barrier()
{
  comm->barrier();
}

template <typename real_t, class shrdmem_class>
void slv_parallel_distmem_mpi<real_t, shrdmem_class>::sndrcv(int peer, int cnt, real_t *ibuf, real_t *obuf)
{
  mpi::request reqs[2];
  int r = 0;
  // isend() / irecv() are non-blocking
  reqs[r++] = comm->isend(peer, msg_halo, obuf, cnt); 
  reqs[r++] = comm->irecv(peer, msg_halo, ibuf, cnt);
  mpi::wait_all(reqs, reqs + r);
}
#endif

// explicit instantiations
#include "cfg/cfg_types.hpp"
#if defined(USE_FLOAT)
template class slv_parallel_distmem_mpi<float, slv_parallel_serial<float>>;
#endif
#if defined(USE_DOUBLE)
template class slv_parallel_distmem_mpi<double, slv_parallel_serial<double>>;
#endif
#if defined(USE_LDOUBLE)
template class slv_parallel_distmem_mpi<long double, slv_parallel_serial<long double>>;
#endif
