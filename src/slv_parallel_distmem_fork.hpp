/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef SLV_PARALLEL_DISTMEM_FORK_HPP
#  define SLV_PARALLEL_DISTMEM_FORK_HPP

#  include "slv_parallel_distmem.hpp"

// fork(), getpid(), pid_t, wait(), _exit()
extern "C" {
#  include <sys/types.h>
#  include <unistd.h>
#  include <sys/wait.h>
};

template <typename real_t, class shrdmem_class>
class slv_parallel_distmem_fork : public slv_parallel_distmem<real_t, shrdmem_class>
{
  private: int fork_init(int nk)
  {
    for (int k = 0; k < nk; ++k)
    {
      // fork
      if (k != 0)
      {
        pid_t pid = fork();
        if (pid == -1) error_macro("fork() failed!")
        else if (pid == 0) return k; // child process
      }
    }
    return 0;
  }

  private: int size, rank;

  public: slv_parallel_distmem_fork(adv<real_t> *fllbck, adv<real_t> *advsch, 
    out<real_t> *output, vel<real_t> *velocity, ini<real_t> *intcond,
    int nx, int ny, int nz, grd<real_t> *grid, quantity<si::time, real_t> dt, int nsd
  ) 
    : slv_parallel_distmem<real_t, shrdmem_class>(
      fllbck, advsch, output, velocity, intcond, nx, ny, nz, grid, dt, size=nsd, rank=fork_init(nsd)
    )
  { }

  public: ~slv_parallel_distmem_fork()
  {
    if (rank == 0)
    {
      int status;
      for (int k = 1; k < size; ++k) if (wait(&status) < 0) error_macro("wait() failed!");
    }
    else _exit(EXIT_SUCCESS);
  }

  private: void distmem_barrier()
  {
cerr << rank << "->barrier()" << endl;
  }

  private: void sndrcv(int peer, int cnt, real_t *ibuf, real_t *obuf)
  {
    /*
    mpi::request reqs[2];
    int r = 0;
    // isend() / irecv() are non-blocking
    reqs[r++] = comm->isend(peer, msg_halo, obuf, cnt); 
    reqs[r++] = comm->irecv(peer, msg_halo, ibuf, cnt);
    mpi::wait_all(reqs, reqs + r);
    */
  }
};

#endif
