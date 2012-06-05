/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef OPT_GRD_HPP
#  define OPT_GRD_HPP

#  include "opt.hpp"
#  include "grd_carthesian.hpp"

inline void opt_grd_desc(po::options_description &desc)
{
  desc.add_options()
    ("grd", po::value<string>()->default_value("carthesian"), "grid: carthesian")
    ("grd.dx", po::value<string>(), "gridbox length (X) [m]")
    ("grd.dy", po::value<string>()->default_value("1"), "gridbox length (Y) [m]")
    ("grd.dz", po::value<string>()->default_value("1"), "gridbox length (Z) [m]")
    ("grd.nx", po::value<int>()->default_value(1), "number of grid points (X)")
    ("grd.ny", po::value<int>()->default_value(1), "number of grid points (Y)")
    ("grd.nz", po::value<int>()->default_value(1), "number of grid points (Z)");
}

template <typename real_t>
grd<real_t> *opt_grd(const po::variables_map& vm)
{
  string grdtype = vm.count("grd") ? vm["grd"].as<string>() : "<unspecified>";

  if (grdtype == "carthesian")
  {
    if (!vm.count("grd.dx")) error_macro("specifying grd.dx is mandatory")
    int 
      nx = vm["grd.nx"].as<int>(),
      ny = vm["grd.ny"].as<int>(),
      nz = vm["grd.nz"].as<int>();

    quantity<si::length, real_t> 
      dx = real_cast<real_t>(vm, "grd.dx") * si::metres,
      dy = real_cast<real_t>(vm, "grd.dy") * si::metres,
      dz = real_cast<real_t>(vm, "grd.dz") * si::metres;

    // some sanity checks
    if (dx / si::metres <= 0 || dy / si::metres <= 0 || dz / si::metres <= 0)
      error_macro("grd.dx, grd.dy, grd.dz must all be > 0") 
    if (nx <= 0 || ny <= 0 || nz <= 0) 
      error_macro("grd.nx, grd.ny and grd.nz must all be > 0")

    return new grd_carthesian<real_t>(dx, dy, dz, nx, ny, nz);
  }
  else error_macro("unsupported grid type: " << grdtype)
}

#endif
