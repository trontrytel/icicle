/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef SLV_VEL_HPP
#  define SLV_VEL_HPP

#  include "opt.hpp"
#  include "slv_parallel_openmp.hpp"
#  include "slv_parallel_threads.hpp"
#  include "slv_parallel_serial.hpp"
#  include "slv_parallel_distmem_mpi.hpp"
#  include "slv_parallel_distmem_fork.hpp"

template <typename real_t>
slv<real_t> *opt_slv(const po::variables_map& vm,
  adv<real_t> *fllbck, 
  adv<real_t> *advsch,
  out<real_t> *output,
  vel<real_t> *velocity,
  ini<real_t> *intcond,
  int nx, int ny, int nz, 
  grd<real_t> *grid,
  quantity<si::time, real_t> dt
)
{
  string slvtype = vm.count("slv") ? vm["slv"].as<string>() : "<unspecified>";
  if (slvtype == "serial")
    return new slv_serial<real_t>(fllbck, advsch, output, velocity, intcond,
      0, nx - 1, nx,
      0, ny - 1, ny,
      0, nz - 1, nz,
      grid, dt
    );
  else
  {
    if (!vm.count("nsd")) error_macro("subdomain count not specified (--nsd option)")
    int nsd = vm["nsd"].as<int>();
#  ifdef _OPENMP
    if (slvtype == "openmp")
      return new slv_parallel_openmp<real_t>(
        fllbck, advsch, output, velocity, intcond, 
        0, nx - 1, nx, 
        0, ny - 1, ny,
        0, nz - 1, nz, 
        grid, dt, nsd
      );
    else
#  endif
#  ifdef USE_BOOST_THREAD
    if (slvtype == "threads")
      return new slv_parallel_threads<real_t>(
        fllbck, advsch, output, velocity, intcond, 
        0, nx - 1, nx, 
        0, ny - 1, ny,
        0, nz - 1, nz, 
        grid, dt, nsd
      );
    else
#  endif
    if (slvtype == "fork")
      return new slv_parallel_distmem_fork<real_t, slv_parallel_serial<real_t> >(
        fllbck, advsch, output, velocity, intcond, nx, ny, nz, grid, dt, nsd
      );
    else
    // TODO: fork+threads, fork+openmp
#  ifdef USE_BOOST_MPI
    if (slvtype == "mpi")
      return new slv_parallel_distmem_mpi<real_t, slv_parallel_serial<real_t> >(
        fllbck, advsch, output, velocity, intcond, nx, ny, nz, grid, dt, nsd
      );
    else
#    ifdef USE_BOOST_THREAD
    if (slvtype == "mpi+threads")
      error_macro("TODO: mpi+threads not implemented yet...")
    else
#    endif
#    ifdef USE_OPENMP
    if (slvtype == "mpi+opemp")
      error_macro("TODO: mpi+openmp not implemented yet...")
    else
#    endif
#  endif
    error_macro("unsupported solver type: " << slvtype)
  }
}

#endif
