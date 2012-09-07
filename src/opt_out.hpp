/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#pragma once

#  include "opt.hpp"
#  include "grd.hpp"
#  include "stp.hpp"
#  include "out_debug.hpp"
#  include "out_gnuplot.hpp"
#  include "out_netcdf.hpp"

inline void opt_out_desc(po::options_description &desc)
{
  desc.add_options()
    ("out", po::value<string>(), "output: debug, gnuplot, netcdf")
    ("out.gnuplot.using", po::value<string>()->default_value("0:-2:1"), "using column specification: 0:-2:1 or 1:-2:2")
    ("out.netcdf.file", po::value<string>()->default_value(""), "output filename")
    ("out.netcdf.ver", po::value<int>()->default_value(4), "netCDF format version (3 or 4)");
}

template <typename real_t>
out<real_t> *opt_out(
  const po::variables_map &vm, 
  const stp<real_t> &setup,  
  const string &cmdline)
{
  string outtype = vm.count("out") ? vm["out"].as<string>() : "<unspecified>";

  // gnuplot
  if (outtype == "gnuplot")
    return new out_gnuplot<real_t>(setup.grid, vm["out.gnuplot.using"].as<string>());
  else

  // debug
  if (outtype == "debug")
    return new out_debug<real_t>();
  else

  // netcdf
  if (outtype == "netcdf")
    return new out_netcdf<real_t>(vm["out.netcdf.file"].as<string>(), setup, 
      vm["out.netcdf.ver"].as<int>(), cmdline);
  else

  error_macro("unsupported output type: " << outtype)
}
