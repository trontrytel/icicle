/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @section DESCRIPTION
 *    a fork()/IPC implementation of slv_parallel_distmem built for several reasons:
 *    - as a facility for testing the distmem logic when compiling without MPI
 *    - for validating that the distmem logic is free of MPI-related assumptions
 *    - for fun :)
 */
#pragma once
#include "slv_parallel_distmem.hpp"

#include <boost/interprocess/ipc/message_queue.hpp>
namespace ipc = boost::interprocess;

template <typename real_t, class shrdmem_class>
class slv_parallel_distmem_fork : public slv_parallel_distmem<real_t, shrdmem_class>
{
  private: string *ipckey;

  private: int fork_init(int nk);

  private: int size, rank;
  private: unique_ptr<ipc::message_queue> *queues_real;
  private: unique_ptr<ipc::message_queue> *queues_bool;

  private: string queue_name(string pfx, int k);
 
  public: slv_parallel_distmem_fork(
    const stp<real_t> &setup, 
    out<real_t> &output, 
    int nsd
  ); 

  public: ~slv_parallel_distmem_fork();

  private: void distmem_barrier();

  private: int qid(int from, int to);

  private: void sndrcv(int peer, int cnt, real_t *ibuf, real_t *obuf);
};
