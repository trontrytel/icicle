/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date January - March 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef OPT_EQS_HPP
#  define OPT_EQS_HPP

#  include "opt.hpp"
#  include "eqs_scalar_advection.hpp"
#  include "eqs_shallow_water.hpp"
#  include "eqs_isentropic.hpp"
#  include "eqs_harmonic_oscillator.hpp"
#  include "eqs_todo_bulk.hpp"

inline void opt_eqs_desc(po::options_description &desc)
{
  desc.add_options()
    ("eqs", po::value<string>()->default_value("scalar_advection"), "equation system: shallow_water, isentropic, ...")

    ("eqs.isentropic.nlev", po::value<int>(), "number of fluid layers")
    ("eqs.isentropic.p_top", po::value<string>()->default_value("0"), "pressure at the uppermost surface [Pa]")
    ("eqs.isentropic.theta_frst", po::value<string>(), "mid-layer potential temperature of the first layer [K]")
    ("eqs.isentropic.abslev", po::value<int>(), "absorber lowermost level")
    ("eqs.isentropic.absamp", po::value<string>()->default_value("1"), "absorber amplitude [1]") 

    ("eqs.harmonic_oscillator.omega", po::value<string>(), "omega [Hz]") ;
  // TODO: constants container --cst.g --cst.cp ...
}

template <typename real_t>
eqs<real_t> *opt_eqs(const po::variables_map& vm, const grd<real_t> &grid, const ini<real_t> &intcond)
{
  string initype= vm.count("eqs") ? vm["eqs"].as<string>() : "<unspecified>";
#  if defined(USE_EQS_SCALAR_ADVECTION)
  if (initype == "scalar_advection")
    return new eqs_scalar_advection<real_t>();
  else 
#  endif
#  if defined(USE_EQS_SHALLOW_WATER)
  if (initype == "shallow_water")
    return new eqs_shallow_water<real_t>(grid);
  else
#  endif
#  if defined(USE_EQS_ISENTROPIC)
  if (initype == "isentropic")
  {
    if (!vm.count("eqs.isentropic.nlev")) error_macro("TODO")
    if (!vm.count("eqs.isentropic.abslev")) error_macro("TODO")
    return new eqs_isentropic<real_t>(grid, 
      vm["eqs.isentropic.nlev"].as<int>(),
      real_cast<real_t>(vm, "eqs.isentropic.p_top") * si::pascals,
      real_cast<real_t>(vm, "eqs.isentropic.theta_frst") * si::kelvins,
      vm["eqs.isentropic.abslev"].as<int>(),
      real_cast<real_t>(vm, "eqs.isentropic.absamp")
    );
  }
  else
#  endif
#  if defined(USE_EQS_HARMONIC_OSCILLATOR)
  if (initype == "harmonic_oscillator")
    return new eqs_harmonic_oscillator<real_t>(
      real_cast<real_t>(vm, "eqs.harmonic_oscillator.omega") / si::seconds
    );
  else 
#  endif
#  if defined(USE_EQS_TODO)
  if (initype == "todo_bulk")
    return new eqs_todo_bulk<real_t>(grid);
  else 
#  endif
  error_macro("unsupported equation system: " << initype)
}

#endif
