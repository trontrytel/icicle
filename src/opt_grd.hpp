/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef OPT_GRD_HPP
#  define OPT_GRD_HPP

#  include "opt.hpp"
#  include "grd_arakawa-c-lorenz.hpp"

#  ifdef ICICLE_OPT_DESCS 
void opt_grd_desc(po::options_description &desc)
{
  desc.add_options()
    ("grd", po::value<string>()->default_value("arakawa-c-lorenz"), "grid: arakawa-c-lorenz")
    ("grd.dx", po::value<string>(), "gridbox length (X) [m]")
    ("grd.dy", po::value<string>()->default_value("1"), "gridbox length (Y) [m]")
    ("grd.dz", po::value<string>()->default_value("1"), "gridbox length (Z) [m]")
    ("grd.nx", po::value<int>()->default_value(1), "number of grid points (X)")
    ("grd.ny", po::value<int>()->default_value(1), "number of grid points (Y)")
    ("grd.nz", po::value<int>()->default_value(1), "number of grid points (Z)");
}
#  endif

template <typename real_t>
grd<real_t> *opt_grd(const po::variables_map& vm)
{
  string grdtype = vm.count("grd") ? vm["grd"].as<string>() : "<unspecified>";

  if (grdtype == "arakawa-c-lorenz")
  {
    if (!vm.count("grd.dx")) error_macro("grd.dx, grd.dy, grd.dz options are mandatory")
    int 
      nx = vm["grd.nx"].as<int>(),
      ny = vm["grd.ny"].as<int>(),
      nz = vm["grd.nz"].as<int>();

    quantity<si::length, real_t> 
      dx = boost::lexical_cast<real_t>(vm["grd.dx"].as<string>()) * si::metres,
      dy = boost::lexical_cast<real_t>(vm["grd.dy"].as<string>()) * si::metres,
      dz = boost::lexical_cast<real_t>(vm["grd.dz"].as<string>()) * si::metres;

    // some sanity checks
    if (dx / si::metres <= 0 || dy / si::metres <= 0 || dz / si::metres <= 0)
      error_macro("grd.dx, grd.dy, grd.dz must all be > 0") 
    if (nx <= 0 || ny <= 0 || nz <= 0) 
      error_macro("grd.nx, grd.ny and grd.nz must all be > 0")

    return new grd_arakawa_c_lorenz<real_t>(dx, dy, dz, nx, ny, nz);
  }
  else error_macro("unsupported grid type: " << grdtype)
}

#endif
