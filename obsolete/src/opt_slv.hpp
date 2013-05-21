/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#pragma once

#  include "opt.hpp"
#  include "slv_parallel_serial.hpp"
#  include "slv_parallel_distmem_mpi.hpp"
#  include "slv_parallel_distmem_fork.hpp"

inline void opt_slv_desc(po::options_description &desc)
{
  desc.add_options()
    ("slv", po::value<string>(), "solver: serial ...")
    ("nsd", po::value<int>(), "subdomain number [1]"); // TODO: rename it to slv.nsd
}

template <typename real_t>
slv<real_t> *opt_slv(
  const po::variables_map &vm, 
  const stp<real_t> &setup, 
  out<real_t> &output
)
{
  string slvtype = vm.count("slv") ? vm["slv"].as<string>() : "<unspecified>";
  if (slvtype == "serial")
    return new slv_parallel_serial<real_t>(setup, output,
      0, setup.grid.nx() - 1, 
      0, setup.grid.ny() - 1, 
      0, setup.grid.nz() - 1
    );
  else if (slvtype != "<unspecified>")
  {
    if (!vm.count("nsd")) error_macro("subdomain count not specified (--nsd option)")
    int nsd = vm["nsd"].as<int>();
    if (slvtype == "fork")
      return new slv_parallel_distmem_fork<real_t, slv_parallel_serial<real_t>>(
        setup, output, nsd
      );
    else
    // TODO: fork+threads, fork+openmp
#  ifdef USE_BOOST_MPI
    if (slvtype == "mpi")
      return new slv_parallel_distmem_mpi<real_t, slv_parallel_serial<real_t>>(
        setup, output, nsd
      );
    else
#  endif
    error_macro("unsupported solver type: " << slvtype)
  }
  else error_macro("unsupported solver type: " << slvtype)
}