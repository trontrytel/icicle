/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date March 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef EQS_TODO_HPP
#  define EQS_TODO_HPP

#  include "cmn.hpp"
#  include "eqs.hpp"
#  include "rhs_explicit.hpp"
#  include "grd.hpp"
#  include "phc.hpp"

template <typename real_t>
class eqs_todo : public eqs<real_t> 
{
  // nested class (well... struct)
  protected: struct params
  {
    quantity<si::mass_density, real_t> rho_unit;
    int idx_rhod, idx_rhod_rv, idx_rhod_th;
  };

  // ctor
  public: eqs_todo(const grd<real_t> &grid, params *par)
  {
    par->rho_unit = 1 * si::kilograms / si::cubic_metres;

    // auxiliary variable
    this->aux.push_back(new struct eqs<real_t>::axv({
      "rhod", "dry air density", this->quan2str(par->rho_unit),
      vector<int>({0, 0, 0})
    }));
    par->idx_rhod = this->aux.size() - 1;

    // eg. eqn 1b in Szumowski, Grabowski & Ochs 1998, Atmos. Res.
    // cf. eqn 3.55 in the Jacobson's Fundamentals of Atmos. Modelling (2nd edn, 2005)
    this->sys.push_back(new struct eqs<real_t>::gte({
      "rhod_rv", "dry air density times water vapour mixing ratio (i.e. water vapour mass density or absolute humidity)",
      this->quan2str(par->rho_unit), 
      typename eqs<real_t>::positive_definite(true),
    }));
    par->idx_rhod_rv = this->sys.size() - 1;

    // eg. eqn 1a in Szumowski, Grabowski & Ochs 1998, Atmos. Res.
    // cf. eqn 3.68 in the Jacobson's Fundamentals of Atmos. Modelling (2nd edn, 2005)
    this->sys.push_back(new struct eqs<real_t>::gte({
      "rhod_th", "dry air density times potential temperature (i.e. energy mass density over specific heat capacity)",
      this->quan2str(par->rho_unit * si::kelvins), 
      typename eqs<real_t>::positive_definite(true),
    }));
    par->idx_rhod_th = this->sys.size() - 1;
  }
};
#endif
