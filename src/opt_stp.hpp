/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date January 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef OPT_STP_HPP
#  define OPT_STP_HPP

#  include "opt.hpp"
#  include "stp.hpp"

inline void opt_stp_desc(po::options_description &desc)
{
  desc.add_options()
    ("dt_out", po::value<string>(), "output interval [s]")
    ("t_max", po::value<string>(), "time of simulation [s]")
    ("dt", po::value<string>()->default_value("auto"), "time step [s]");
}

template <typename real_t>
stp<real_t> *opt_stp(const po::variables_map& vm, 
  adv<real_t> *advsch,
  adv<real_t> *fllbck,
  vel<real_t> *velocity,
  ini<real_t> *intcond,
  grd<real_t> *grid,
  eqs<real_t> *equations
  )
{ 
  quantity<si::time, real_t> dt=0*si::seconds;
  if (vm["dt"].as<string>() != "auto" ) dt = real_cast<real_t>(vm, "dt") * si::seconds; 
  return new stp<real_t>(
    fllbck, advsch,
    velocity,
    intcond,
    grid,
    equations,
    real_cast<real_t>(vm, "dt_out")*si::seconds,
    real_cast<real_t>(vm, "t_max")*si::seconds,
    dt
  );
}

#endif
