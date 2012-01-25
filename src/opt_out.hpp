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
#  include "stp.hpp"
#  include "eqs.hpp"
#  include "out_debug.hpp"
#  include "out_gnuplot.hpp"
#  include "out_netcdf.hpp"

#  ifdef ICICLE_OPT_DESCS 
void opt_out_desc(po::options_description &desc)
{
  desc.add_options()
    ("out", po::value<string>(), "output: debug, gnuplot, netcdf")
    ("out.gnuplot.using", po::value<string>()->default_value("0:-2:1"), "using column specification: 0:-2:1 or 1:-2:2")
    ("out.netcdf.file", po::value<string>(), "output filename (e.g. for netcdf)")
    ("out.netcdf.ver", po::value<int>()->default_value(4), "netCDF format version (3 or 4)");
}
#  endif

template <typename real_t>
out<real_t> *opt_out(const po::variables_map &vm, 
  stp<real_t> *setup, grd<real_t> *grid, eqs<real_t> &equation,
  const string &options)
{
  string outtype = vm.count("out") ? vm["out"].as<string>() : "<unspecified>";
  if (outtype == "gnuplot")
    return new out_gnuplot<real_t>(grid, vm["out.gnuplot.using"].as<string>());
  else
  if (outtype == "debug")
    return new out_debug<real_t>();
  else
#  ifdef USE_NETCDF
  if (outtype == "netcdf")
  {
    if (!vm.count("out.netcdf.file")) error_macro("output filename not specified (--out.netcdf.file option)")
    return new out_netcdf<real_t>(vm["out.netcdf.file"].as<string>(), setup, grid, equation, 
      vm["out.netcdf.ver"].as<int>(), options);
  }
  else
#  endif
  error_macro("unsupported output type: " << outtype)
}

#endif
