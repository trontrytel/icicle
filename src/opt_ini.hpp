/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef OPT_INI_HPP
#  define OPT_INI_HPP

#  include "opt.hpp"
#  include "ini_origin-one.hpp"
#  include "ini_cone.hpp"

template <typename real_t>
ini<real_t> *opt_ini(const po::variables_map& vm, grd<real_t> *grid)
{
  string initype= vm.count("ini") ? vm["ini"].as<string>() : "<unspecified>";
  if (initype == "origin-one")
  {
    return new ini_origin_one<real_t>(grid->dx(), grid->dy(), grid->dz());
  }
  else if (initype == "cone")
  {
    if (!vm.count("ini.cone.h") || !vm.count("ini.cone.x0") || !vm.count("ini.cone.z0") || !vm.count("ini.cone.r"))
      error_macro("ini.cone.[h,x0,z0,r] must be specified")
    quantity<si::length, real_t> 
      h  = boost::lexical_cast<real_t>(vm["ini.cone.h" ].as<string>()) * si::metres,
      x0 = boost::lexical_cast<real_t>(vm["ini.cone.x0"].as<string>()) * si::metres,
      z0 = boost::lexical_cast<real_t>(vm["ini.cone.z0"].as<string>()) * si::metres,
      r  = boost::lexical_cast<real_t>(vm["ini.cone.r" ].as<string>()) * si::metres;
    return new ini_cone<real_t>(h, x0, z0, r);
  }
  else error_macro("unsupported initial condition: " << initype)
}

#endif
