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
#  include "slv_serial.hpp"
#  include "slv_parallel_openmp.hpp"
#  include "slv_parallel_threads.hpp"
#  include "ini.hpp" // TODO TEMP

template <typename real_t>
slv<si::dimensionless, real_t> *opt_slv(const po::variables_map& vm,
  adv<si::dimensionless, real_t> *fllbck, 
  adv<si::dimensionless, real_t> *advsch,
  out<si::dimensionless, real_t> *output,
  vel<real_t> *velocity,
  ini<real_t> *intcond,
  int nx, int ny, int nz, 
  grd<real_t> *grid,
  quantity<si::time, real_t> dt
)
{
  string slvtype = vm.count("slv") ? vm["slv"].as<string>() : "<unspecified>";
  if (slvtype == "serial")
    return new slv_serial<si::dimensionless, real_t>(fllbck, advsch, output, velocity, intcond,
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
      return new slv_parallel_openmp<si::dimensionless, real_t>(
        fllbck, advsch, output, velocity, intcond, nx, ny, nz, grid, dt, nsd);
    else
#  endif
#  ifdef USE_BOOST_THREAD
    if (slvtype == "threads")
      return new slv_parallel_threads<si::dimensionless, real_t>(
        fllbck, advsch, output, velocity, intcond, nx, ny, nz, grid, dt, nsd);
    else
#  endif
    error_macro("unsupported solver type: " << slvtype)
  }
}

#endif
