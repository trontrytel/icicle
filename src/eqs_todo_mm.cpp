/** @file
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date August 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 @  @brief contains definitions of eqs_todo_mm class - 2 moment parametrisation of warm rain  microphysics
 */

#include "eqs_todo_mm.hpp"
#include "phc.hpp"
#include "phc_theta.hpp"
#include "phc_kelvin_term.hpp"
#include "phc_mm.hpp"

/** @brief the 2 moment parametrisation of warm rain microphysics
 *
 * Consult Morrison and Grabowski 2007.
 * 
 * predicted variables: concentrations and mixing ratios of cloud and drizzle/rain water
 */

template <typename real_t>
template <bool fill_tmp>

class eqs_todo_mm<real_t>::activation : public rhs_explicit<real_t>
{
  private: struct params &par;
  private: real_t mean_rd;
  private: real_t sdev_rd;
  private: real_t n_tot;
  private: real_t chem_b;
  private: real_t scl;

  // ctor
  public: activation(
    struct params &par,
    real_t scl,
    real_t mean_rd, 
    real_t sdev_rd,
    real_t n_tot,
    real_t chem_b
  )
  : par(par), scl(scl), mean_rd(mean_rd), sdev_rd(sdev_rd), n_tot(n_tot), chem_b(chem_b)
  {}

  public: void explicit_part(
    mtx::arr<real_t> &Ra,
    ptr_unordered_map<string, mtx::arr<real_t>> &aux,
    const mtx::arr<real_t> * const * const psi,
    const quantity<si::time, real_t> dt
  ) const
  {
    mtx::arr<real_t>
      &tmp = aux.at("tmp");
    const mtx::arr<real_t>
      &rhod = aux.at("rhod"),
      &rhod_rv = (*psi[par.idx_rhod_rv]),     
      &rhod_th = (*psi[par.idx_rhod_th]),        
      &rhod_rl = (*psi[par.idx_rhod_rl]),
      &rhod_nl = (*psi[par.idx_rhod_nl]);

    if (fill_tmp)
    {
      for (int k = rhod.lbound(mtx::k); k <= rhod.ubound(mtx::k); ++k)
      {
        for (int i = rhod.lbound(mtx::i); i <= rhod.ubound(mtx::i); ++i)
        {
          for (int j = rhod.lbound(mtx::j); j <= rhod.ubound(mtx::j); ++j)
          {
          //TODO beta as an option (instead of 0.5) ?
          tmp(i,j,k) = rhod(i,j,k) * std::max<real_t>(0, (   
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
          }
        }
      }        
    }
    Ra(Ra.ijk) += scl * tmp(Ra.ijk); 
  }
};

template <typename real_t>
class eqs_todo_mm<real_t>::cond_evap : public rhs_explicit<real_t>
{
  private: struct params &par;
  private: real_t scl;

  // ctor
  public: cond_evap(struct params &par, real_t scl)
  : par(par), scl(scl)
  {}

  public: void explicit_part(
    mtx::arr<real_t> &Rce,
    ptr_unordered_map<string, mtx::arr<real_t>> &aux,
    const mtx::arr<real_t> * const * const psi,
    const quantity<si::time, real_t> dt
  ) const
  {
    mtx::arr<real_t>
      &tmp = aux.at("tmp");
    const mtx::arr<real_t>
      &rhod = aux.at("rhod"),
      &rhod_rv = (*psi[par.idx_rhod_rv]),     
      &rhod_th = (*psi[par.idx_rhod_th]),        
      &rhod_rl = (*psi[par.idx_rhod_rl]),
      &rhod_nl = (*psi[par.idx_rhod_nl]),
      &rhod_rr = (*psi[par.idx_rhod_rr]),
      &rhod_nr = (*psi[par.idx_rhod_nr]);

//TODO rr and nr for rain cond/evap
      for (int k = rhod.lbound(mtx::k); k <= rhod.ubound(mtx::k); ++k)
      {
        for (int i = rhod.lbound(mtx::i); i <= rhod.ubound(mtx::i); ++i)
        {
          for (int j = rhod.lbound(mtx::j); j <= rhod.ubound(mtx::j); ++j)
          {

          Rce(i,j,k) += 
          rhod(i,j,k) * ((rhod_rv(i,j,k) / rhod(i,j,k) 
          -
          phc::r_vs<real_t>(
            phc::T<real_t>(
              rhod_th(i,j,k) / rhod(i,j,k) * si::kelvins, 
              phc::p<real_t>(rhod_th(i,j,k) * si::kelvins * si::kilograms / si::cubic_metres, rhod_rv(i,j,k) / rhod(i,j,k)), 
              rhod_rv(i,j,k) / rhod(i,j,k)), 
            phc::p<real_t>(rhod_th(i,j,k) * si::kelvins * si::kilograms / si::cubic_metres, rhod_rv(i,j,k) / rhod(i,j,k))
          )
          ) / 
          (phc::mm::tau_relax<real_t>( 
            real_t(pow(rhod_rl(i,j,k) / rhod_nl(i,j,k) / 
              (real_t(4./3 * M_PI) * (phc::rho_w<real_t>() / si::kilograms *si::cubic_metres)), real_t(1./3))) * si::metres,  
              //TODO why do I need real_t(pow)*si::   ? 
            rhod_nl(i,j,k) / rhod(i,j,k) / si::cubic_metres) 
           / si::seconds)
          /
          (
          real_t(1) + 
          real_t(1.) / (rhod(i,j,k) * si::kilograms / si::cubic_metres) / phc::R_v<real_t>() / 
            phc::T<real_t>(
              rhod_th(i,j,k) / rhod(i,j,k) * si::kelvins, 
              phc::p<real_t>(rhod_th(i,j,k) * si::kelvins * si::kilograms / si::cubic_metres, rhod_rv(i,j,k) / rhod(i,j,k)), 
              rhod_rv(i,j,k) / rhod(i,j,k)
            ) * 
          phc::mm::desdT<real_t>(
            phc::p_vs<real_t>(
              phc::T<real_t>(
                rhod_th(i,j,k) / rhod(i,j,k) * si::kelvins, 
                phc::p<real_t>(rhod_th(i,j,k) * si::kelvins * si::kilograms / si::cubic_metres, rhod_rv(i,j,k) / rhod(i,j,k)), 
                rhod_rv(i,j,k) / rhod(i,j,k))), 
            phc::T<real_t>(
              rhod_th(i,j,k) / rhod(i,j,k) * si::kelvins, 
              phc::p<real_t>(rhod_th(i,j,k) * si::kelvins * si::kilograms / si::cubic_metres, rhod_rv(i,j,k) / rhod(i,j,k)), 
              rhod_rv(i,j,k) / rhod(i,j,k))
          ) * 
          phc::l_v<real_t>(
            phc::T<real_t>(
             rhod_th(i,j,k) / rhod(i,j,k) * si::kelvins, 
              phc::p<real_t>(rhod_th(i,j,k) * si::kelvins * si::kilograms / si::cubic_metres, rhod_rv(i,j,k) / rhod(i,j,k)), 
              rhod_rv(i,j,k) / rhod(i,j,k))
          ) / phc::c_p<real_t>(rhod_rv(i,j,k) / rhod(i,j,k))
          ));
          }
        }
      }
  }
};

// ctor
template <typename real_t>
eqs_todo_mm<real_t>::eqs_todo_mm(
  const grd<real_t> &grid, 
  map<enum processes, bool> opts, 
  real_t mean_rd, 
  real_t sdev_rd,
  real_t n_tot,
  real_t chem_b
)
  : eqs_todo<real_t>(grid, &par), opts(opts)
{ 
  // in activation term all activated droples are assumed to have the radius of 1 um
  const real_t ccnmass = real_t(4./3 * M_PI * 1e-18) * phc::rho_w<real_t>() / si::kilograms * si::cubic_metres;

  if (opts[act]) this->sys.at(par.idx_rhod_rv).rhs_terms.push_back(new activation<true>(par, -ccnmass, mean_rd, sdev_rd, n_tot, chem_b));
  if (opts[cond]) this->sys.at(par.idx_rhod_rv).rhs_terms.push_back(new cond_evap(par, 0));

  this->sys.push_back(new struct eqs<real_t>::gte({
    "rhod_rl", "dry air density times liquid water mixing ratio (i.e. liquid water mass density)",
    this->quan2str(this->par.rho_unit),
    eqs<real_t>::positive_definite(true),
  }));
  if (opts[act]) this->sys.back().rhs_terms.push_back(new activation<false>(par, ccnmass, mean_rd, sdev_rd, n_tot, chem_b));
  par.idx_rhod_rl = this->sys.size() - 1;

  this->sys.push_back(new struct eqs<real_t>::gte({
    "rhod_nl", "dry air density times liquid water number concentration",
    this->quan2str(this->par.rho_unit),
    eqs<real_t>::positive_definite(true),
  }));
  if (opts[act]) this->sys.back().rhs_terms.push_back(new activation<false>(par, 1, mean_rd, sdev_rd, n_tot, chem_b));
  par.idx_rhod_nl = this->sys.size() - 1;

  this->sys.push_back(new struct eqs<real_t>::gte({
    "rhod_rr", "dry air density times rain water mixing ratio (i.e. rain water mass density)",
    this->quan2str(this->par.rho_unit),
    eqs<real_t>::positive_definite(true),
  }));
  if (opts[cond]) this->sys.back().rhs_terms.push_back(new cond_evap(par, 0));
  par.idx_rhod_rr = this->sys.size() - 1;

  this->sys.push_back(new struct eqs<real_t>::gte({
    "rhod_nr", "dry air density times rain water number concentration",
    this->quan2str(this->par.rho_unit),
    eqs<real_t>::positive_definite(true),
  }));
  if (opts[cond]) this->sys.back().rhs_terms.push_back(new cond_evap(par, 0));
  par.idx_rhod_nr = this->sys.size() - 1;

  // auxliary variable for temporary data
  ptr_map_insert(this->aux)("tmp", typename eqs<real_t>::axv({
    "tmp", "temporary data", "n/a",
    typename eqs<real_t>::invariable(false),
    vector<int>({0, 0, 0})
  })); // TODO: flag it as not to be output
}

// explicit instantiations
#include "cfg/cfg_types.hpp"
#if defined(USE_FLOAT)
template class eqs_todo_mm<float>;
#endif
#if defined(USE_DOUBLE)
template class eqs_todo_mm<double>;
#endif
#if defined(USE_LDOUBLE)
template class eqs_todo_mm<long double>;
#endif
