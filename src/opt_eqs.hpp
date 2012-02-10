/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date January 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef OPT_EQS_HPP
#  define OPT_EQS_HPP

#  include "opt.hpp"
#  include "eqs_scalar_advection.hpp"
#  include "eqs_shallow_water_2d.hpp"

inline void opt_eqs_desc(po::options_description &desc)
{
  desc.add_options()
    ("eqs", po::value<string>()->default_value("scalar_advection"), "equation system: shallow_water_2d");
}

template <typename real_t>
eqs<real_t> *opt_eqs(const po::variables_map& vm)
{
  string initype= vm.count("eqs") ? vm["eqs"].as<string>() : "<unspecified>";
  if (initype == "scalar_advection")
    return new eqs_scalar_advection<real_t>();
  else if (initype == "shallow_water_2d")
    return new eqs_shallow_water_2d<real_t>();
  else error_macro("unsupported equation system: " << initype)
}

#endif
