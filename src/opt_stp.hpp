/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date January 2012
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef OPT_STP_HPP
#  define OPT_STP_HPP

#  include "opt.hpp"
#  include "stp.hpp"

#  ifdef ICICLE_OPT_DESCS 
void opt_stp_desc(po::options_description &desc)
{
  desc.add_options()
    ("dt_out", po::value<string>(), "output interval [s]")
    ("t_max", po::value<string>(), "time of simulation [s]")
    ("dt", po::value<string>()->default_value("auto"), "time step [s]");
}
#  endif

template <typename real_t>
stp<real_t> *opt_stp(const po::variables_map& vm, 
  adv<real_t> *advsch,
  adv<real_t> *fllbck,
  vel<real_t> *velocity,
  ini<real_t> *intcond,
  grd<real_t> *grid,
  int nx, int ny, int nz
  )
{ 
  quantity<si::time, real_t> dt=0*si::seconds;
  if (vm["dt"].as<string>() != "auto" ) dt = boost::lexical_cast<real_t>(vm["dt"].as<string>()) * si::seconds; 
  if (!vm.count("t_max")) error_macro("t_max not specified")
  if (!vm.count("dt_out")) error_macro("dt_out not specified")
  return new stp<real_t>(
    fllbck, advsch,
    velocity,
    intcond,
    grid,
    nx, ny, nz,
    boost::lexical_cast<real_t>(vm["dt_out"].as<string>())*si::seconds,
    boost::lexical_cast<real_t>(vm["t_max"].as<string>())*si::seconds,
    dt
  );
}

#endif
