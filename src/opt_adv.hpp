/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef OPT_ADV_HPP
#  define OPT_ADV_HPP

#  include "opt.hpp"
#  include "adv_mpdata.hpp"
#  include "adv_mpdata_fct.hpp"
#  include "adv_leapfrog.hpp"
#  include "adv_lax-wendroff.hpp"

template <typename real_t>
void opt_adv(const po::variables_map& vm,
  adv<real_t> **fllbck,
  adv<real_t> **advsch, 
  grd<real_t> *grid
)
{
  grd_arakawa_c_lorenz<real_t> *g = dynamic_cast<grd_arakawa_c_lorenz<real_t>*>(grid);
  if (g == NULL) error_macro("all advection schemes are compatible with the Arakawa-C grid only!")

  string advscheme = vm.count("adv") ? vm["adv"].as<string>() : "<unspecified>";

  *fllbck = NULL;
  if (advscheme == "leapfrog")
  {
    *advsch = new adv_leapfrog<real_t>(g);
    *fllbck = new adv_mpdata<real_t>(g, 1);
  }
  else if (advscheme == "mpdata")
  {
    if (vm["adv.mpdata.fct"].as<bool>())
      *advsch = new adv_mpdata_fct<real_t>(g, vm["adv.mpdata.iord"].as<int>());
    else
      *advsch = new adv_mpdata<real_t>(g, vm["adv.mpdata.iord"].as<int>());
  }
  else error_macro("unsupported advection scheme: " << advscheme)
}

#endif
