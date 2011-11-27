/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef OPT_OUT_HPP
#  define OPT_OUT_HPP

#  include "opt.hpp"
#  include "grd.hpp"
#  include "out_gnuplot.hpp"
#  include "out_netcdf.hpp"

template <typename real_t>
out<si::dimensionless, real_t> *opt_out(const po::variables_map& vm, 
  grd<real_t> *grid, int nx, int ny, int nz)
{
  string outtype = vm.count("out") ? vm["out"].as<string>() : "<unspecified>";
  if (outtype == "gnuplot")
    return new out_gnuplot<si::dimensionless, real_t>();
  else
#  ifdef USE_NETCDF
  if (outtype == "netcdf")
  {
    if (!vm.count("out.netcdf.file")) error_macro("output filename not specified (--out.netcdf.file option)")
    return new out_netcdf<si::dimensionless, real_t>(vm["out.netcdf.file"].as<string>(), grid, nx, ny, nz);
  }
  else
#  endif
  error_macro("unsupported output type: " << outtype)
}

#endif
