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
#  include "eqs_shallow_water.hpp"
#  include "eqs_isentropic.hpp"

inline void opt_eqs_desc(po::options_description &desc)
{
  desc.add_options()
    ("eqs", po::value<string>()->default_value("scalar_advection"), "equation system: shallow_water, isentropic, ...")
    ("eqs.isentropic.nlev", po::value<int>(), "number of fluid layers")
    ("eqs.isentropic.g", po::value<string>()->default_value("9.81"), "acceleration due to gravity [m/s2]");
  // TODO: constants container --cst.g --cst.cp ...
}

template <typename real_t>
eqs<real_t> *opt_eqs(const po::variables_map& vm, const grd<real_t> &grid, const ini<real_t> &intcond)
{
  string initype= vm.count("eqs") ? vm["eqs"].as<string>() : "<unspecified>";
  if (initype == "scalar_advection")
    return new eqs_scalar_advection<real_t>();
  else if (initype == "shallow_water")
    return new eqs_shallow_water<real_t>(grid, intcond);
  else if (initype == "isentropic")
  {
    if (!vm.count("eqs.isentropic.nlev")) error_macro("TODO")
    return new eqs_isentropic<real_t>(grid, intcond,
      vm["eqs.isentropic.nlev"].as<int>(),
      real_cast<real_t>(vm, "eqs.isentropic.g") * si::metres_per_second_squared
    );
  }
  else error_macro("unsupported equation system: " << initype)
}

#endif
