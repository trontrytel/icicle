/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#include "slv_parallel_serial.hpp"
#include "slv_parallel_distmem_fork.hpp"

// fork(), getpid(), pid_t, wait(), _exit()
extern "C" {
#  include <sys/types.h>
#  include <unistd.h>
#  include <sys/wait.h>
};

template <typename real_t, class shrdmem_class>
int slv_parallel_distmem_fork<real_t, shrdmem_class>::fork_init(int nk)
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

template <typename real_t, class shrdmem_class>
string slv_parallel_distmem_fork<real_t, shrdmem_class>::queue_name(string pfx, int k)
{
    ostringstream tmp;
    tmp << *ipckey << "." << pfx << "." << k; 
    return tmp.str();
}
 
template <typename real_t, class shrdmem_class>
slv_parallel_distmem_fork<real_t, shrdmem_class>::slv_parallel_distmem_fork(
    const stp<real_t> &setup, 
    out<real_t> &output, 
    int nsd
) 
: slv_parallel_distmem<real_t, shrdmem_class>(setup, output, size=nsd, rank=fork_init(nsd))
{ 
    queues_real = new unique_ptr<ipc::message_queue>[2 * nsd];
    queues_bool = new unique_ptr<ipc::message_queue>[nsd];
    for (int kk = 0; kk < 2 * size; ++kk) 
    {   
      queues_real[kk].reset(new ipc::message_queue(ipc::open_or_create, queue_name("real", kk).c_str(), 
        1, (setup.advsch.stencil_extent() - 1) / 2 * setup.grid.nz() * setup.grid.ny() * sizeof(real_t))
      ); 
    }
    for (int kk = 0; kk < size; ++kk) 
    {
      queues_bool[kk].reset(new ipc::message_queue(ipc::open_or_create, queue_name("bool", kk).c_str(), 
        kk == 0 ? size - 1 : 1, sizeof(bool))
      ); 
    }  
}

template <typename real_t, class shrdmem_class>
slv_parallel_distmem_fork<real_t, shrdmem_class>::~slv_parallel_distmem_fork()
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

template <typename real_t, class shrdmem_class>
void slv_parallel_distmem_fork<real_t, shrdmem_class>::distmem_barrier()
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

template <typename real_t, class shrdmem_class>
int slv_parallel_distmem_fork<real_t, shrdmem_class>::qid(int from, int to)
{
    assert(((to + 1 + size) % size == from) || ((to - 1 + size) % size == from));
    return to + (((to + 1) % size == from) ? 1 : 0);
}

template <typename real_t, class shrdmem_class>
void slv_parallel_distmem_fork<real_t, shrdmem_class>::sndrcv(int peer, int cnt, real_t *ibuf, real_t *obuf)
{
    cnt *= sizeof(real_t);
    queues_real[qid(rank, peer)]->send(obuf, cnt, 0);
    size_t rcvd;
    unsigned int prio;
    queues_real[qid(peer, rank)]->receive(ibuf, cnt, rcvd, prio);
    assert(cnt == rcvd && prio == 0);
}

// explicit instantiations
#include "cfg/cfg_types.hpp"
#if defined(USE_FLOAT)
template class slv_parallel_distmem_fork<float, slv_parallel_serial<float>>;
#endif
#if defined(USE_DOUBLE)
template class slv_parallel_distmem_fork<double, slv_parallel_serial<double>>;
#endif
#if defined(USE_LDOUBLE)
template class slv_parallel_distmem_fork<long double, slv_parallel_serial<long double>>;
#endif
