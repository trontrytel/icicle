/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#pragma once

#  include "opt.hpp"
#  include "adv_upstream.hpp"
#  include "adv_mpdata.hpp"
#  include "adv_leapfrog.hpp"

inline void opt_adv_desc(po::options_description &desc)
{
  desc.add_options()
    ("adv", po::value<string>(), "advection scheme: leapfrog, upstream, mpdata")
    ("adv.mpdata.iord", po::value<int>()->default_value(2), "MPDATA iteration count: 1, 2, ...")
    ("adv.mpdata.cross_terms", po::value<bool>()->default_value(true), "MPDATA cross terms: 0 (off) or 1 (on)")
    ;
}

template <typename real_t>
adv<real_t> *opt_adv(const po::variables_map& vm)
{
  string advscheme = vm.count("adv") ? vm["adv"].as<string>() : "<unspecified>";

  if (advscheme == "leapfrog")
    return new adv_leapfrog<real_t>();
  else if (advscheme == "upstream")
    return new adv_upstream<real_t>();
  else if (advscheme == "mpdata")
  {
      return new adv_mpdata<real_t>(
        vm["adv.mpdata.iord"].as<int>(),
        vm["adv.mpdata.cross_terms"].as<bool>()
      );
  }
  else error_macro("unsupported advection scheme: " << advscheme)
}
