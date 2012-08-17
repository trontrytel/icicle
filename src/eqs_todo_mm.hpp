/** @file
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date August 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#pragma once
#include "eqs_todo.hpp"
#include "phc.hpp"
#include "phc_theta.hpp"
#include "phc_kelvin_term.hpp"

template <typename real_t>
class eqs_todo_mm : public eqs_todo<real_t> 
{
  // nested class (well... struct)
  protected: struct params : eqs_todo<real_t>::params
  {
    int idx_rhod_rl, idx_rhod_rr, idx_rhod_nl, idx_rhod_nr; // advected variables indices
  };

  protected: params par;

  // a container for storing options (i.e. which processes ar on/off)
  public: enum processes {act, cond, acc, autoc, self, turb, sedi};
  private: map<enum processes, bool> opts;

  // nested class
  private: class activation : public rhs_explicit<real_t>
  {
    private: struct params &par;
    private: real_t sign;
    private: real_t mean_rd;
    private: real_t sdev_rd;
    private: real_t n_tot;
    private: real_t chem_b;

    //ctor
    public: activation(
      struct params &par,
      real_t sign,
      real_t mean_rd, 
      real_t sdev_rd,
      real_t n_tot,
      real_t chem_b
    )
    : par(par), sign(sign), mean_rd(mean_rd), sdev_rd(sdev_rd), n_tot(n_tot), chem_b(chem_b)
    {}

    public: void explicit_part(
      mtx::arr<real_t> &Ra,
      const ptr_unordered_map<string, mtx::arr<real_t>> &aux,
      const mtx::arr<real_t> * const * const psi,
      const quantity<si::time, real_t> dt
    ) const
    {
      const mtx::arr<real_t>
        &rhod = aux.at("rhod"),
        &rhod_rv = (*psi[par.idx_rhod_rv]),     
        &rhod_th = (*psi[par.idx_rhod_th]),        
        &rhod_rl = (*psi[par.idx_rhod_rl]),
        &rhod_nl = (*psi[par.idx_rhod_nl]);

      for (int k = rhod.lbound(mtx::k); k <= rhod.ubound(mtx::k); ++k)
      {
        for (int i = rhod.lbound(mtx::i); i <= rhod.ubound(mtx::i); ++i)
        {
          for (int j = rhod.lbound(mtx::j); j <= rhod.ubound(mtx::j); ++j)
          {
          //TODO beta as an option (instead of 0.5) ?
          Ra(i,j,k) += rhod(i,j,k) * sign * std::max<real_t>(0, (   
            (n_tot / real_t(2.) * erfc(log(
             (pow(mean_rd, -(1+0.5)) * sqrt (real_t(4) * 
              pow(phc::kelvin::A<real_t>(
                phc::T<real_t>(
                  rhod_th(i,j,k) / rhod(i,j,k) * si::kelvins, 
                  phc::p<real_t>(rhod_th(i,j,k) * si::kelvins * si::kilograms / si::cubic_metres, rhod_rv(i,j,k) / rhod(i,j,k)), 
                  rhod_rv(i,j,k) / rhod(i,j,k) 
                )
              ) / si::metres, 3) / real_t(27) / chem_b))
             /
             (rhod_rv(i,j,k) / rhod(i,j,k) / 
              phc::r_vs<real_t>(
                phc::T<real_t>(
                  rhod_th(i,j,k) / rhod(i,j,k) * si::kelvins, 
                  phc::p<real_t>(rhod_th(i,j,k) * si::kelvins * si::kilograms / si::cubic_metres, rhod_rv(i,j,k) / rhod(i,j,k)), 
                  rhod_rv(i,j,k) / rhod(i,j,k)
                ), 
                phc::p<real_t>(rhod_th(i,j,k) * si::kelvins * si::kilograms / si::cubic_metres, rhod_rv(i,j,k) / rhod(i,j,k))
              ) - real_t(1))
             )) 
             / sqrt(real_t(2)) / log(pow(sdev_rd, 1+0.5)))
            -
            rhod_nl(i,j,k) / rhod(i,j,k)) / (dt / si::seconds));

//TODO for rl  and rv should be multimplied by that
//* real_t(4/3) * M_PI * (1e-18) * phc::rho_w<real_t>() / si::kilograms * si::metres

          }
        }
      }        

    }
  };

  // ctor
  public: eqs_todo_mm(
    const grd<real_t> &grid, 
    map<enum processes, bool> opts, 
    real_t mean_rd, 
    real_t sdev_rd,
    real_t n_tot,
    real_t chem_b
  )
    : eqs_todo<real_t>(grid, &par), opts(opts)
  {
    if (opts[act]) this->sys.at(par.idx_rhod_rv).rhs_terms.push_back(new activation(par,-1 ,mean_rd, sdev_rd, n_tot, chem_b));

    this->sys.push_back(new struct eqs<real_t>::gte({
      "rhod_rl", "dry air density times liquid water mixing ratio (i.e. liquid water mass density)",
      this->quan2str(this->par.rho_unit),
      eqs<real_t>::positive_definite(true),
    }));
    if (opts[act]) this->sys.back().rhs_terms.push_back(new activation(par, 1, mean_rd, sdev_rd, n_tot, chem_b));
    par.idx_rhod_rl = this->sys.size() - 1;

    this->sys.push_back(new struct eqs<real_t>::gte({
      "rhod_nl", "dry air density times liquid water number concentration",
      this->quan2str(this->par.rho_unit),
      eqs<real_t>::positive_definite(true),
    }));
    if (opts[act]) this->sys.back().rhs_terms.push_back(new activation(par, 1, mean_rd, sdev_rd, n_tot, chem_b));
    par.idx_rhod_nl = this->sys.size() - 1;
  }
};
