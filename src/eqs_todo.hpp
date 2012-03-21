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

template <typename real_t>
class eqs_todo : public eqs<real_t> 
{
  // nested class (well... struct)
  private: struct params
  {
    quantity<si::mass_density, real_t> rho_unit;
  };

  // provate field
  private: params par;

  // ctor
  public: eqs_todo(const grd<real_t> &grid)
  {
    par.rho_unit = 1 * si::kilograms / si::cubic_metres;

    // eg. eqn 1b in Szumowski, Grabowski & Ochs 1998, Atmos. Res.
    // cf. eqn 3.55 in the Jacobson's Fundamentals of Atmos. Modelling (2nd edn, 2005)
    this->sys.push_back(new struct eqs<real_t>::gte({
      "rho_v", "water vapour mass density (i.e. absolute humidity or dry air density times water vapour mixing ratio)",
      this->quan2str(par.rho_unit), 
      typename eqs<real_t>::positive_definite(true),
    }));
    // cf. eqn 3.68 in the Jacobson's Fundamentals of Atmos. Modelling (2nd edn, 2005)
    this->sys.push_back(new struct eqs<real_t>::gte({
      "E", "energy mass density (i.e. dry air heat capacity times virtual potential temperature times dry air density)",
      this->quan2str(par.rho_unit), 
      typename eqs<real_t>::positive_definite(true),
    }));
    this->init_maps();
  }
};
#endif
