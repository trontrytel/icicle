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

template <typename real_t>
grd<real_t> *opt_grd(const po::variables_map& vm)
{
  string grdtype = vm.count("grd") ? vm["grd"].as<string>() : "<unspecified>";

  if (grdtype == "arakawa-c")
  {
    if (!vm.count("grd.dx") || !vm.count("grd.dy") || !vm.count("grd.dz"))
      error_macro("grd.dx, grd.dy, grd.dz options are mandatory")

    quantity<si::length, real_t> 
      dx = boost::lexical_cast<real_t>(vm["grd.dx"].as<string>()) * si::metres,
      dy = boost::lexical_cast<real_t>(vm["grd.dy"].as<string>()) * si::metres,
      dz = boost::lexical_cast<real_t>(vm["grd.dz"].as<string>()) * si::metres;

    // some sanity checks
   if (dx / si::metres <= 0 || dy / si::metres <= 0 || dz / si::metres <= 0)
     error_macro("grd.dx, grd.dy, grd.dz must all be >= 0") 
   assert(dx == dy && dy == dz); // TODO...

    return new grd_arakawa_c_lorenz<real_t>(dx, dy, dz);
  }
  else error_macro("unsupported grid type: " << grdtype)
}

#endif
