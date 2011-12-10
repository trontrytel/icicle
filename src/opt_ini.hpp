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
#  include "grd.hpp"
#  include "ini_boxcar.hpp"
#  include "ini_pks_wwg_1989.hpp"
#  include "ini_cone.hpp"

void opt_ini_desc(po::options_description &desc)
{
  desc.add_options()
    ("ini", po::value<string>(), "initical condition: boxcar, cone, pks_wwg_1989")
    ("ini.boxcar.a", po::value<string>()->default_value("0"), "a [m]")
    ("ini.boxcar.b", po::value<string>(), "b [m]")
    ("ini.boxcar.A", po::value<string>()->default_value("1"), "A [1]")
    ("ini.boxcar.A0", po::value<string>()->default_value("0"), "A0 [1]")
    ("ini.cone.h", po::value<string>(), "h [m]")
    ("ini.cone.x0", po::value<string>(), "x0 [m]")
    ("ini.cone.z0", po::value<string>(), "z0 [m]")
    ("ini.cone.r", po::value<string>(), "r [m]");
}

template <typename real_t>
ini<real_t> *opt_ini(const po::variables_map& vm, grd<real_t> *grid)
{
  string initype= vm.count("ini") ? vm["ini"].as<string>() : "<unspecified>";
  if (initype == "pks_wwg_1989")
    return new ini_pks_wwg_1989<real_t>(grid->dx());
  else if (initype == "boxcar")
  {
    quantity<si::length, real_t> 
      a = boost::lexical_cast<real_t>(vm["ini.boxcar.a"].as<string>()) * si::metres,
      b = boost::lexical_cast<real_t>(vm["ini.boxcar.b"].as<string>()) * si::metres;
    quantity<si::dimensionless, real_t>
      A  = boost::lexical_cast<real_t>(vm["ini.boxcar.A"].as<string>()),
      A0 = boost::lexical_cast<real_t>(vm["ini.boxcar.A0"].as<string>());
    return new ini_boxcar<real_t>(a, b, A, A0);
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
