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
  private: struct params
  {
    quantity<si::mass_density, real_t> rho_unit;
    int idx_rhod, idx_rhod_rv, idx_rhod_th, idx_rhod_rl; // auxiliary variable indices
  };

  // private field
  private: params par;

  // TODO: should it be virtual? and place in another file???
  // the saturation adjustment (aka ,,bulk'' microphysics)
  public: void adjustments(
    int n,
    const ptr_vector<mtx::arr<real_t>> &aux, 
    const vector<ptr_vector<mtx::arr<real_t>>> &psi
  ) 
  {
    const mtx::arr<real_t> 
      &rhod = aux[par.idx_rhod],
      &rhod_rv = psi[par.idx_rhod_rv][n],
      &rhod_rl = psi[par.idx_rhod_rl][n],
      &rhod_th = psi[par.idx_rhod_th][n];

    for (int i = rhod.lbound(mtx::i); i <= rhod.ubound(mtx::i); ++i)
      for (int j = rhod.lbound(mtx::j); j <= rhod.ubound(mtx::j); ++j)
        for (int k = rhod.lbound(mtx::k); k <= rhod.ubound(mtx::k); ++k)
        {
          quantity<phc::mixing_ratio, real_t> 
            r = rhod_rv(i,j,k) / rhod(i,j,k);
          quantity<si::pressure, real_t> 
            p = phc::p_1000<real_t>() * real_t(pow(
              (rhod_th(i,j,k) * phc::R_d<real_t>() * si::kelvins * si::kilograms / si::cubic_metres) 
                / phc::p_1000<real_t>() * (real_t(1) + r / phc::eps<real_t>()),
              1 / (1 - phc::R_over_c_p(r))
            ));
          quantity<si::temperature, real_t> 
            T = rhod_th(i,j,k) / rhod(i,j,k) * si::kelvins * phc::exner<real_t>(p, r);

          cerr << "RH(" << i << "," << j << "," << k << ")=" 
            << phc::r_vs<real_t>(T,p) / r << endl;
        }
  }

  // ctor
  public: eqs_todo(const grd<real_t> &grid)
  {
    par.rho_unit = 1 * si::kilograms / si::cubic_metres;

    // auxiliary variable
    this->aux.push_back(new struct eqs<real_t>::axv({
      "rhod", "dry air density", this->quan2str(par.rho_unit),
      vector<int>({0, 0, 0})
    }));
    par.idx_rhod = this->aux.size() - 1;

    // eg. eqn 1b in Szumowski, Grabowski & Ochs 1998, Atmos. Res.
    // cf. eqn 3.55 in the Jacobson's Fundamentals of Atmos. Modelling (2nd edn, 2005)
    this->sys.push_back(new struct eqs<real_t>::gte({
      "rhod_rv", "dry air density times water vapour mixing ratio (i.e. water vapour mass density or absolute humidity)",
      this->quan2str(par.rho_unit), 
      typename eqs<real_t>::positive_definite(true),
    }));
    par.idx_rhod_rv = this->sys.size() - 1;

    // only for bulk model (i.e. not for super droplets)
    if (true) // TODO!
    {
      this->sys.push_back(new struct eqs<real_t>::gte({
        "rhod_rl", "dry air density times liquid water mixing ratio (i.e. liquid water mass density)",
        this->quan2str(par.rho_unit), 
        typename eqs<real_t>::positive_definite(true),
      }));
      par.idx_rhod_rl = this->sys.size() - 1;
    }

    // eg. eqn 1a in Szumowski, Grabowski & Ochs 1998, Atmos. Res.
    // cf. eqn 3.68 in the Jacobson's Fundamentals of Atmos. Modelling (2nd edn, 2005)
    this->sys.push_back(new struct eqs<real_t>::gte({
      "rhod_th", "dry air density times potential temperature (i.e. energy mass density over specific heat capacity)",
      this->quan2str(par.rho_unit * si::kelvins), 
      typename eqs<real_t>::positive_definite(true),
    }));
    par.idx_rhod_th = this->sys.size() - 1;
  }
};
#endif
