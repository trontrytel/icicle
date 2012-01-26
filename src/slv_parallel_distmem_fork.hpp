/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 *  @section DESCRIPTION
 *    a fork()/IPC implementation of slv_parallel_distmem built for several reasons:
 *    - as a facility for testing the distmem logic when compiling without MPI
 *    - for validating that the distmem logic is free of MPI-related assumptions
 *    - for fun :)
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

#  include <boost/interprocess/ipc/message_queue.hpp>
namespace ipc = boost::interprocess;

template <typename real_t, class shrdmem_class>
class slv_parallel_distmem_fork : public slv_parallel_distmem<real_t, shrdmem_class>
{
  private: string *ipckey;

  private: int fork_init(int nk)
  {
    for (int k = 0; k < nk; ++k)
    {
      // a unique ipc key - must be set before fork()
      if (k == 0)  
      {   
        ostringstream tmp;
        tmp << "icicle." << getpid(); 
        ipckey = new string(tmp.str()); 
      } 
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
  private: auto_ptr<ipc::message_queue> *queues_real;
  private: auto_ptr<ipc::message_queue> *queues_bool;

  private: string queue_name(string pfx, int k)
  {
    ostringstream tmp;
    tmp << *ipckey << "." << pfx << "." << k; 
    return tmp.str();
  }
 
  public: slv_parallel_distmem_fork(stp<real_t> *setup, out<real_t> *output, int nsd
  ) 
    : slv_parallel_distmem<real_t, shrdmem_class>(setup, output, size=nsd, rank=fork_init(nsd))
  { 
    queues_real = new auto_ptr<ipc::message_queue>[2 * nsd];
    queues_bool = new auto_ptr<ipc::message_queue>[nsd];
    for (int kk = 0; kk < 2 * size; ++kk) 
    {   
      queues_real[kk].reset(new ipc::message_queue(ipc::open_or_create, queue_name("real", kk).c_str(), 
        1, (setup->advsch->stencil_extent() - 1) / 2 * setup->grid->nz() * setup->grid->ny() * sizeof(real_t))
      ); 
    }
    for (int kk = 0; kk < size; ++kk) 
    {
      queues_bool[kk].reset(new ipc::message_queue(ipc::open_or_create, queue_name("bool", kk).c_str(), 
        kk == 0 ? size - 1 : 1, sizeof(bool))
      ); 
    }  
  }

  public: ~slv_parallel_distmem_fork()
  {
    delete[] queues_real;
    delete[] queues_bool;
    if (rank == 0)
    {
      for (int kk = 0; kk < 2 * size; ++kk) ipc::message_queue::remove(queue_name("real", kk).c_str());
      for (int kk = 0; kk < size; ++kk) ipc::message_queue::remove(queue_name("bool", kk).c_str());
    }
    delete ipckey;
    if (rank == 0)
    {
      int status;
      for (int k = 1; k < size; ++k) if (wait(&status) < 0) error_macro("wait() failed!");
    }
    else _exit(EXIT_SUCCESS);
  }

  private: void distmem_barrier()
  {
    bool buf;
    size_t rcvd;
    unsigned int prio;
    if (rank == 0)
    {   
      for (int k=1; k < size; ++k) queues_bool[0]->receive(&buf, sizeof(bool), rcvd, prio);
      for (int k=1; k < size; ++k) queues_bool[k]->send(&buf, sizeof(bool), 0);
    }   
    else 
    {   
      queues_bool[0]->send(&buf, sizeof(bool), 0);
      queues_bool[rank]->receive(&buf, sizeof(bool), rcvd, prio);
    }   
  }

  private: int qid(int from, int to)
  {
    assert(((to + 1 + size) % size == from) || ((to - 1 + size) % size == from));
    return to + (((to + 1) % size == from) ? 1 : 0);
  }

  private: void sndrcv(int peer, int cnt, real_t *ibuf, real_t *obuf)
  {
    cnt *= sizeof(real_t);
    queues_real[qid(rank, peer)]->send(obuf, cnt, 0);
    size_t rcvd;
    unsigned int prio;
    queues_real[qid(peer, rank)]->receive(ibuf, cnt, rcvd, prio);
    assert(cnt == rcvd && prio == 0);
  }
};

#endif
