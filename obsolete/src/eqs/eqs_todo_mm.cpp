/** @file
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date October 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 @  @brief contains definitions of eqs_todo_mm class - 2 moment parametrisation of warm rain  microphysics
 */

#include "eqs_todo_mm.hpp"
#include "../phc/phc_kessler.hpp"
#include "../phc/phc_theta.hpp"
#include "../phc/phc_kelvin_term.hpp"
#include "../phc/phc_mm.hpp"

/** @brief the 2 moment parametrisation of warm rain microphysics
 * Consult Morrison and Grabowski 2007.
 * predicted variables: concentrations and mixing ratios of cloud and drizzle/rain water
 */

template <typename real_t>
template <bool fill_tmp_act>
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
      &tmp_act = aux.at("tmp_act");
    const mtx::arr<real_t>
      &rhod = aux.at("rhod"),
      &rhod_rv = (*psi[par.idx_rhod_rv]),     
      &rhod_th = (*psi[par.idx_rhod_th]),        
      &rhod_rl = (*psi[par.idx_rhod_rl]),
      &rhod_nl = (*psi[par.idx_rhod_nl]);

    if (fill_tmp_act)
    {
      for (int k = rhod.lbound(mtx::k); k <= rhod.ubound(mtx::k); ++k)
      {
        for (int j = rhod.lbound(mtx::j); j <= rhod.ubound(mtx::j); ++j)
        {
          for (int i = rhod.lbound(mtx::i); i <= rhod.ubound(mtx::i); ++i)
          {
            //TODO beta as an option (instead of 0.5) ?
  
            //helper temperature and pressure fields
            quantity<si::pressure, real_t> p = phc::p<real_t>(
              rhod_th(i,j,k) * si::kelvins * si::kilograms / si::cubic_metres, 
              rhod_rv(i,j,k) / rhod(i,j,k));
            quantity<si::temperature, real_t> T = phc::T<real_t>(
              rhod_th(i,j,k) / rhod(i,j,k) * si::kelvins, p, rhod_rv(i,j,k) / rhod(i,j,k));

            tmp_act(i,j,k) = rhod(i,j,k) * std::max<real_t>(0, (   
              (n_tot / real_t(2.) * 
                erfc( std::log( (std::pow( mean_rd, -(1+0.5)) * sqrt (real_t(4) * 
                  std::pow(phc::kelvin::A<real_t>(T) / si::metres, 3) / real_t(27) / chem_b)) /
                  (rhod_rv(i,j,k) / rhod(i,j,k) / phc::r_vs<real_t>(T, p) - real_t(1)))) /  
                sqrt(real_t(2)) / std::log(std::pow(sdev_rd, 1+0.5))) -
              rhod_nl(i,j,k) / rhod(i,j,k)) / (dt / si::seconds));
          }
        }
      }        
    }
    Ra(Ra.ijk) += scl * tmp_act(Ra.ijk); 
  }
};

template <typename real_t>
class eqs_todo_mm<real_t>::activation_th : public rhs_explicit<real_t>
{
  private: struct params &par;
  private: real_t scl;

  // ctor
  public: activation_th(
    struct params &par,
    real_t scl
  )
  : par(par), scl(scl)
  {}

  public: void explicit_part(
    mtx::arr<real_t> &Ra_th,
    ptr_unordered_map<string, mtx::arr<real_t>> &aux,
    const mtx::arr<real_t> * const * const psi,
    const quantity<si::time, real_t> 
  ) const
  {
    mtx::arr<real_t>
      &tmp_act = aux.at("tmp_act");
    const mtx::arr<real_t>
      &rhod = aux.at("rhod"),
      &rhod_rv = (*psi[par.idx_rhod_rv]),     
      &rhod_th = (*psi[par.idx_rhod_th]);        

    for (int k = rhod.lbound(mtx::k); k <= rhod.ubound(mtx::k); ++k)
    {
      for (int j = rhod.lbound(mtx::j); j <= rhod.ubound(mtx::j); ++j)
      {
        for (int i = rhod.lbound(mtx::i); i <= rhod.ubound(mtx::i); ++i)
        {
          //TODO beta as an option (instead of 0.5) ?

          //helper temperature and pressure fields
          quantity<si::pressure, real_t> p = phc::p<real_t>(
            rhod_th(i,j,k) * si::kelvins * si::kilograms / si::cubic_metres, 
            rhod_rv(i,j,k) / rhod(i,j,k));
          quantity<si::temperature, real_t> T = phc::T<real_t>(
            rhod_th(i,j,k) / rhod(i,j,k) * si::kelvins, p, rhod_rv(i,j,k) / rhod(i,j,k));

          Ra_th(i,j,k) +=  scl * phc::dtheta_drv<real_t>(T, p, rhod_rv(i,j,k) / rhod(i,j,k), 
              rhod_th(i,j,k) * si::kelvins * si::kilograms / si::cubic_metres, 
              rhod(i,j,k) * si::kilograms / si::cubic_metres) / si::kelvins *
            tmp_act(i,j,k); 
        }
      }
    }        
  }
};

template <typename real_t>
template <bool fill_tmp_cond>
class eqs_todo_mm<real_t>::cond_evap : public rhs_explicit<real_t>
{
  private: struct params &par;
  private: real_t sgn;

  // ctor
  public: cond_evap(struct params &par, real_t sgn)
  : par(par), sgn(sgn)
  {}

  public: void explicit_part(
    mtx::arr<real_t> &Rce,
    ptr_unordered_map<string, mtx::arr<real_t>> &aux,
    const mtx::arr<real_t> * const * const psi,
    const quantity<si::time, real_t> dt
  ) const
  {
    mtx::arr<real_t>
      &tmp_cond = aux.at("tmp_cond");
    const mtx::arr<real_t>
      &rhod = aux.at("rhod"),
      &rhod_rv = (*psi[par.idx_rhod_rv]),     
      &rhod_th = (*psi[par.idx_rhod_th]),        
      &rhod_rl = (*psi[par.idx_rhod_rl]),
      &rhod_nl = (*psi[par.idx_rhod_nl]);

//TODO if conditions for cloud/no cloud, rain/no rain conditions
//TODO rr and nr for rain cond/evap
//TODO change in the number of cloud and rain droplets due to evaporation

    for (int k = rhod.lbound(mtx::k); k <= rhod.ubound(mtx::k); ++k)
    {
      for (int j = rhod.lbound(mtx::j); j <= rhod.ubound(mtx::j); ++j)
      {
        for (int i = rhod.lbound(mtx::i); i <= rhod.ubound(mtx::i); ++i)
        {
          if (fill_tmp_cond)
          {
            // helper temperature and pressure fields
            quantity<si::pressure, real_t> p = phc::p<real_t>(
              rhod_th(i,j,k) * si::kelvins * si::kilograms / si::cubic_metres, 
              rhod_rv(i,j,k) / rhod(i,j,k));
            quantity<si::temperature, real_t> T = phc::T<real_t>(
              rhod_th(i,j,k) / rhod(i,j,k) * si::kelvins, p, rhod_rv(i,j,k) / rhod(i,j,k));
  
            quantity<si::dimensionless, real_t> helper;
  
            if (rhod_rl(i,j,k) == 0); // do nothing
            else
            { 
              tmp_cond(i,j,k) = rhod(i,j,k) * (rhod_rv(i,j,k) / rhod(i,j,k) - phc::r_vs<real_t>(T,p)) / 
                (phc::mm::tau_relax<real_t>( 
                  std::pow(rhod_rl(i,j,k) / rhod_nl(i,j,k) / 
                    (real_t(4./3 * M_PI) * (phc::rho_w<real_t>() / si::kilograms *si::cubic_metres)), real_t(1./3)) * si::metres,  
                  rhod_nl(i,j,k) / rhod(i,j,k) / si::cubic_metres) 
                 / si::seconds) / 
               (real_t(1.) + 
                  real_t(1.) / (rhod(i,j,k) * si::kilograms / si::cubic_metres) / phc::R_v<real_t>() / T *
                  phc::mm::desdT<real_t>(phc::p_vs<real_t>(T),T) * 
                  phc::l_v<real_t>(T) / phc::c_p<real_t>(rhod_rv(i,j,k) / rhod(i,j,k))
               );
            }
          }
          if (real_t(-1) * tmp_cond(i,j,k) > rhod_rl(i,j,k)) 
            Rce(i,j,k) += - real_t(.5) * rhod_rl(i,j,k) / (dt / si::seconds);
          else 
            Rce(i,j,k) += sgn * tmp_cond(i,j,k);                          
        }           // the only places with rhod_rl are tau_relax and then helper condition to avaoid negative values
      }             // should number of the droplet changes be linked with this helper condition (how it is done in KK2000)?
    }
  }
};


template <typename real_t>
template <bool fill_tmp_ce_rain>
class eqs_todo_mm<real_t>::ce_rain : public rhs_explicit<real_t>
{
  private: struct params &par;
  private: real_t sgn;

  // ctor
  public: ce_rain(struct params &par, real_t sgn)
  : par(par), sgn(sgn)
  {}

  public: void explicit_part(
    mtx::arr<real_t> &Rcer,
    ptr_unordered_map<string, mtx::arr<real_t>> &aux,
    const mtx::arr<real_t> * const * const psi,
    const quantity<si::time, real_t> dt
  ) const
  {
    mtx::arr<real_t>
      &tmp_ce_rain = aux.at("tmp_ce_rain");
    const mtx::arr<real_t>
      &rhod = aux.at("rhod"),
      &rhod_rv = (*psi[par.idx_rhod_rv]),     
      &rhod_th = (*psi[par.idx_rhod_th]),        
      &rhod_rr = (*psi[par.idx_rhod_rr]),
      &rhod_nr = (*psi[par.idx_rhod_nr]);

//TODO if conditions for cloud/no cloud, rain/no rain conditions
//TODO change in the number of cloud and rain droplets due to evaporation

//TODO no ventilation coeffs -> ok for drizzle but not so much for rain
    for (int k = rhod.lbound(mtx::k); k <= rhod.ubound(mtx::k); ++k)
    {
      for (int j = rhod.lbound(mtx::j); j <= rhod.ubound(mtx::j); ++j)
      {
        for (int i = rhod.lbound(mtx::i); i <= rhod.ubound(mtx::i); ++i)
        {
          if (fill_tmp_ce_rain)
          {
            // helper temperature and pressure fields
            quantity<si::pressure, real_t> p = phc::p<real_t>(
              rhod_th(i,j,k) * si::kelvins * si::kilograms / si::cubic_metres, 
              rhod_rv(i,j,k) / rhod(i,j,k));
            quantity<si::temperature, real_t> T = phc::T<real_t>(
              rhod_th(i,j,k) / rhod(i,j,k) * si::kelvins, p, rhod_rv(i,j,k) / rhod(i,j,k));
  
            quantity<si::dimensionless, real_t> helper;
  
            if (rhod_rr(i,j,k) == 0); // do nothing
            else
            { 
              tmp_ce_rain(i,j,k) = rhod(i,j,k) * (rhod_rv(i,j,k) / rhod(i,j,k) - phc::r_vs<real_t>(T,p)) / 
                (phc::mm::tau_relax<real_t>( 
                  std::pow(rhod_rr(i,j,k) / rhod_nr(i,j,k) / 
                    (real_t(4./3 * M_PI) * (phc::rho_w<real_t>() / si::kilograms *si::cubic_metres)), real_t(1./3)) * si::metres,  
                  rhod_nr(i,j,k) / rhod(i,j,k) / si::cubic_metres) 
                 / si::seconds) / 
               (real_t(1.) + 
                  real_t(1.) / (rhod(i,j,k) * si::kilograms / si::cubic_metres) / phc::R_v<real_t>() / T *
                  phc::mm::desdT<real_t>(phc::p_vs<real_t>(T),T) * 
                  phc::l_v<real_t>(T) / phc::c_p<real_t>(rhod_rv(i,j,k) / rhod(i,j,k))
               );
            }
          }
          if (real_t(-1) * tmp_ce_rain(i,j,k) > rhod_rr(i,j,k)) 
            Rcer(i,j,k) += - real_t(.5) * rhod_rr(i,j,k) / (dt / si::seconds);
          else 
            Rcer(i,j,k) += sgn * tmp_ce_rain(i,j,k);                          
        }           // the only places with rhod_rr are tau_relax and then helper condition to avaoid negative values
      }             // should number of the droplet changes be linked with this helper condition (how it is done in KK2000)?
    }
  }
};

template <typename real_t>
class eqs_todo_mm<real_t>::n_evap : public rhs_explicit<real_t>
{
  private: struct params &par;
  // ctor
  public: n_evap(struct params &par)
  : par(par)
  {}

  public: void explicit_part(
    mtx::arr<real_t> &Rne,
    ptr_unordered_map<string, mtx::arr<real_t>> &aux,
    const mtx::arr<real_t> * const * const psi,
    const quantity<si::time, real_t> dt
  ) const
  {
    mtx::arr<real_t>
      &tmp_cond = aux.at("tmp_cond");
    const mtx::arr<real_t>
      &rhod = aux.at("rhod"),
      &rhod_rl = (*psi[par.idx_rhod_rl]),
      &rhod_nl = (*psi[par.idx_rhod_nl]);

    for (int k = rhod.lbound(mtx::k); k <= rhod.ubound(mtx::k); ++k)
    {
      for (int j = rhod.lbound(mtx::j); j <= rhod.ubound(mtx::j); ++j)
      {
        for (int i = rhod.lbound(mtx::i); i <= rhod.ubound(mtx::i); ++i)
        {
          if (real_t(-1) * tmp_cond(i,j,k) > rhod_rl(i,j,k)) Rne(i,j,k) += - real_t(.5) * rhod_nl(i,j,k) / (dt / si::seconds); 
        }           
      }            
    }
  }
};

template <typename real_t>
class eqs_todo_mm<real_t>::nr_evap : public rhs_explicit<real_t>
{
  private: struct params &par;
  // ctor
  public: nr_evap(struct params &par)
  : par(par)
  {}

  public: void explicit_part(
    mtx::arr<real_t> &Rnre,
    ptr_unordered_map<string, mtx::arr<real_t>> &aux,
    const mtx::arr<real_t> * const * const psi,
    const quantity<si::time, real_t> dt
  ) const
  {
    mtx::arr<real_t>
      &tmp_ce_rain = aux.at("tmp_ce_rain");
    const mtx::arr<real_t>
      &rhod = aux.at("rhod"),
      &rhod_rr = (*psi[par.idx_rhod_rr]),
      &rhod_nr = (*psi[par.idx_rhod_nr]);

    for (int k = rhod.lbound(mtx::k); k <= rhod.ubound(mtx::k); ++k)
    {
      for (int j = rhod.lbound(mtx::j); j <= rhod.ubound(mtx::j); ++j)
      {
        for (int i = rhod.lbound(mtx::i); i <= rhod.ubound(mtx::i); ++i)
        {
          if (real_t(-1) * tmp_ce_rain(i,j,k) > real_t(0)) Rnre(i,j,k) += 
            rhod(i,j,k) * rhod_nr(i,j,k) / rhod_rr(i,j,k) * tmp_ce_rain(i,j,k); 
        }           
      }            
    }
  }
};


template <typename real_t>
class eqs_todo_mm<real_t>::theta_cond_evap : public rhs_explicit<real_t>
{
  private: struct params &par;
  // ctor
  public: theta_cond_evap(struct params &par)
  : par(par)
  {}

  public: void explicit_part(
    mtx::arr<real_t> &Rtce,
    ptr_unordered_map<string, mtx::arr<real_t>> &aux,
    const mtx::arr<real_t> * const * const psi,
    const quantity<si::time, real_t> dt
  ) const
  {
    mtx::arr<real_t>
      &tmp_cond = aux.at("tmp_cond");
    const mtx::arr<real_t>
      &rhod = aux.at("rhod"),
      &rhod_rv = (*psi[par.idx_rhod_rv]),     
      &rhod_th = (*psi[par.idx_rhod_th]),      
      &rhod_rl = (*psi[par.idx_rhod_rl]),
      &rhod_nl = (*psi[par.idx_rhod_nl]);

    for (int k = rhod.lbound(mtx::k); k <= rhod.ubound(mtx::k); ++k)
    {
      for (int j = rhod.lbound(mtx::j); j <= rhod.ubound(mtx::j); ++j)
      {
        for (int i = rhod.lbound(mtx::i); i <= rhod.ubound(mtx::i); ++i)
        {
          quantity<si::pressure, real_t> p = phc::p<real_t>(
            rhod_th(i,j,k) * si::kelvins * si::kilograms / si::cubic_metres, 
            rhod_rv(i,j,k) / rhod(i,j,k));
          quantity<si::temperature, real_t> T = phc::T<real_t>(
            rhod_th(i,j,k) / rhod(i,j,k) * si::kelvins, p, rhod_rv(i,j,k) / rhod(i,j,k));

          if (real_t(-1) * tmp_cond(i,j,k) > rhod_rl(i,j,k))
          {
            Rtce(i,j,k) += - real_t(.5) * rhod_rl(i,j,k) / (dt / si::seconds) * 
              phc::dtheta_drv<real_t>(T, p, rhod_rv(i,j,k) / rhod(i,j,k), 
               rhod_th(i,j,k) * si::kelvins * si::kilograms / si::cubic_metres, 
               rhod(i,j,k) * si::kilograms / si::cubic_metres) / si::kelvins;
          }
          else
          {
            Rtce(i,j,k) += - tmp_cond(i,j,k) * 
              phc::dtheta_drv<real_t>(T, p, rhod_rv(i,j,k) / rhod(i,j,k), 
               rhod_th(i,j,k) * si::kelvins * si::kilograms / si::cubic_metres, 
               rhod(i,j,k) * si::kilograms / si::cubic_metres) / si::kelvins;
          }
        }           
      }            
    }
  }
};

template <typename real_t>
class eqs_todo_mm<real_t>::theta_ce_rain : public rhs_explicit<real_t>
{
  private: struct params &par;
  // ctor
  public: theta_ce_rain(struct params &par)
  : par(par)
  {}

  public: void explicit_part(
    mtx::arr<real_t> &Rtcer,
    ptr_unordered_map<string, mtx::arr<real_t>> &aux,
    const mtx::arr<real_t> * const * const psi,
    const quantity<si::time, real_t> dt
  ) const
  {
    mtx::arr<real_t>
      &tmp_ce_rain = aux.at("tmp_ce_rain");
    const mtx::arr<real_t>
      &rhod = aux.at("rhod"),
      &rhod_rv = (*psi[par.idx_rhod_rv]),     
      &rhod_th = (*psi[par.idx_rhod_th]),      
      &rhod_rr = (*psi[par.idx_rhod_rr]),
      &rhod_nr = (*psi[par.idx_rhod_nr]);

    for (int k = rhod.lbound(mtx::k); k <= rhod.ubound(mtx::k); ++k)
    {
      for (int j = rhod.lbound(mtx::j); j <= rhod.ubound(mtx::j); ++j)
      {
        for (int i = rhod.lbound(mtx::i); i <= rhod.ubound(mtx::i); ++i)
        {
          quantity<si::pressure, real_t> p = phc::p<real_t>(
            rhod_th(i,j,k) * si::kelvins * si::kilograms / si::cubic_metres, 
            rhod_rv(i,j,k) / rhod(i,j,k));
          quantity<si::temperature, real_t> T = phc::T<real_t>(
            rhod_th(i,j,k) / rhod(i,j,k) * si::kelvins, p, rhod_rv(i,j,k) / rhod(i,j,k));

          if (real_t(-1) * tmp_ce_rain(i,j,k) > rhod_rr(i,j,k))
          {
            Rtcer(i,j,k) += - real_t(.5) * rhod_rr(i,j,k) / (dt / si::seconds) * 
              phc::dtheta_drv<real_t>(T, p, rhod_rv(i,j,k) / rhod(i,j,k), 
               rhod_th(i,j,k) * si::kelvins * si::kilograms / si::cubic_metres, 
               rhod(i,j,k) * si::kilograms / si::cubic_metres) / si::kelvins;
          }
          else
          {
            Rtcer(i,j,k) += - tmp_ce_rain(i,j,k) * 
              phc::dtheta_drv<real_t>(T, p, rhod_rv(i,j,k) / rhod(i,j,k), 
               rhod_th(i,j,k) * si::kelvins * si::kilograms / si::cubic_metres, 
               rhod(i,j,k) * si::kilograms / si::cubic_metres) / si::kelvins;
          }
        }           
      }            
    }
  }
};


template <typename real_t>
template <bool fill_tmp_auto>
class eqs_todo_mm<real_t>::autoconv : public rhs_explicit<real_t>
{
  private: struct params &par;
  private: real_t scl;

  // ctor
  public: autoconv(struct params &par, real_t scl)
  : par(par), scl(scl)
  {}

  public: void explicit_part(
    mtx::arr<real_t> &Rauto,
    ptr_unordered_map<string, mtx::arr<real_t>> &aux,
    const mtx::arr<real_t> * const * const psi,
    const quantity<si::time, real_t> dt
  ) const
  {
    mtx::arr<real_t>
      &tmp_auto = aux.at("tmp_auto");
    const mtx::arr<real_t>
      &rhod = aux.at("rhod"),
      &rhod_rv = (*psi[par.idx_rhod_rv]),     
      &rhod_th = (*psi[par.idx_rhod_th]),        
      &rhod_rl = (*psi[par.idx_rhod_rl]),
      &rhod_nl = (*psi[par.idx_rhod_nl]),
      &rhod_rr = (*psi[par.idx_rhod_rr]),
      &rhod_nr = (*psi[par.idx_rhod_nr]);

    if (fill_tmp_auto)
    {
    // TODO: matrix op + where()
      for (int k = rhod.lbound(mtx::k); k <= rhod.ubound(mtx::k); ++k)
      {
        for (int j = rhod.lbound(mtx::j); j <= rhod.ubound(mtx::j); ++j)
        { 
          for (int i = rhod.lbound(mtx::i); i <= rhod.ubound(mtx::i); ++i)
          {
            if (rhod_nl(i,j,k) < 4e-36 || rhod_rl(i,j,k) < 4e-36) tmp_auto(i,j,k) = real_t(0);
            else  // if something is too small e-179 it becomes negative 
                  // so instead of rhod_nl == 0 we have rhod_nl < mass of electron neutrino
            {
              //as in Khairoutdinov and Kogan 2000
              //but taken from Wood 2005 (table 1) - si units!
              //Wood 2005 uses LWC instead of mixing ratios (hence the lack of rhod^-1.47)
              tmp_auto(i,j,k) = rhod(i,j,k) * real_t(7.42*1e13) * 
                std::pow(rhod_rl(i,j,k) / rhod(i,j,k), 2.47) * std::pow(rhod_nl(i,j,k) / rhod(i,j,k), -1.79); 
            }      
          }
        }
      }
    }
    Rauto(Rauto.ijk) += tmp_auto(tmp_auto.ijk) / scl;
  }
};



template <typename real_t>
template <bool fill_tmp_acc>
class eqs_todo_mm<real_t>::accretion : public rhs_explicit<real_t>
{
  private: struct params &par;
  private: real_t scl;

  // ctor
  public: accretion(struct params &par, real_t scl)
  : par(par), scl(scl)
  {}

  public: void explicit_part(
    mtx::arr<real_t> &Racc,
    ptr_unordered_map<string, mtx::arr<real_t>> &aux,
    const mtx::arr<real_t> * const * const psi,
    const quantity<si::time, real_t> dt
  ) const
  {
    mtx::arr<real_t>
      &tmp_acc = aux.at("tmp_acc");
    const mtx::arr<real_t>
      &rhod = aux.at("rhod"),
      &rhod_rv = (*psi[par.idx_rhod_rv]),     
      &rhod_th = (*psi[par.idx_rhod_th]),        
      &rhod_rl = (*psi[par.idx_rhod_rl]),
      &rhod_nl = (*psi[par.idx_rhod_nl]),
      &rhod_rr = (*psi[par.idx_rhod_rr]),
      &rhod_nr = (*psi[par.idx_rhod_nr]);

    if (fill_tmp_acc)
    {
    // TODO: matrix op + where()
      for (int k = rhod.lbound(mtx::k); k <= rhod.ubound(mtx::k); ++k)
      {
        for (int j = rhod.lbound(mtx::j); j <= rhod.ubound(mtx::j); ++j)
        { 
          for (int i = rhod.lbound(mtx::i); i <= rhod.ubound(mtx::i); ++i)
          {
            if (rhod_nl(i,j,k) < 4e-36 || rhod_rl(i,j,k) < 4e-36) tmp_acc(i,j,k) = real_t(0);
            else  // if something is too small e-179 it becomes negative 
                  // so instead of rhod_nl == 0 we have rhod_nl < mass of electron neutrino
            {
              //as in Khairoutdinov and Kogan 2000
              tmp_acc(i,j,k) = rhod(i,j,k) * real_t(67) * 
                std::pow(rhod_rl(i,j,k) * rhod_rr(i,j,k) / rhod(i,j,k) / rhod(i,j,k), 1.15); 
            }      
          }
        }
      }
    }
    Racc(Racc.ijk) += tmp_acc(tmp_acc.ijk) / scl;
  }
};

template <typename real_t>
class eqs_todo_mm<real_t>::col_sink : public rhs_explicit<real_t>
{
  private: struct params &par;
  private: const grd<real_t> &grid;
  // ctor
  public: col_sink(struct params &par, const grd<real_t> &grid)
  : par(par), grid(grid)
  {}

  public: void explicit_part(
    mtx::arr<real_t> &Rcs,
    ptr_unordered_map<string, mtx::arr<real_t>> &aux,
    const mtx::arr<real_t> * const * const psi,
    const quantity<si::time, real_t> dt
  ) const
  {
    mtx::arr<real_t>
      &tmp_auto = aux.at("tmp_auto"),
      &tmp_acc  = aux.at("tmp_acc");
    const mtx::arr<real_t>
      &rhod = aux.at("rhod"),
      &rhod_rl = (*psi[par.idx_rhod_rl]),
      &rhod_nl = (*psi[par.idx_rhod_nl]);

    for (int k = rhod.lbound(mtx::k); k <= rhod.ubound(mtx::k); ++k)
    {
      for (int j = rhod.lbound(mtx::j); j <= rhod.ubound(mtx::j); ++j)
      {
        for (int i = rhod.lbound(mtx::i); i <= rhod.ubound(mtx::i); ++i)
        {
          if (tmp_auto(i,j,k) + tmp_acc(i,j,k) > 4e-36) Rcs(i,j,k) += -
               (tmp_auto(i,j,k) + tmp_acc(i,j,k)) * 
               rhod_nl(i,j,k) / 
               rhod_rl(i,j,k) /
               rhod(i,j,k) / 
               (grid.dx() * grid.dy() * grid.dz() / si::cubic_metres); 
        }           
      }            
    }
  }
};
/*
template <typename real_t>
class eqs_todo_mm<real_t>::terminal_velocity : public rhs_explicit<real_t>
  {
    private: struct params &par;
    private: const grd<real_t> &grid;

    //ctor
    public: terminal_velocity(struct params &par, const grd<real_t> &grid)
    : par(par), grid(grid)
    {}

    public: void explicit_part(
      mtx::arr<real_t> &Rt,
      ptr_unordered_map<string, mtx::arr<real_t>> &aux,
      const mtx::arr<real_t> * const * const psi,
      const quantity<si::time, real_t>
    ) const
    {
      const mtx::arr<real_t>
        &rhod = aux.at("rhod"),
        &rhod_rr = (*psi[par.idx_rhod_rr]);

      assert(min(rhod_rr) >=0);

      for (int k = rhod.lbound(mtx::k); k <= rhod.ubound(mtx::k); ++k)
      {
        for (int i = rhod.lbound(mtx::i); i <= rhod.ubound(mtx::i); ++i)
        {
          quantity<si::velocity, real_t> v_iph = 0 * si::metres_per_second;
          for (int j = rhod.ubound(mtx::j); j > rhod.lbound(mtx::j); --j)
          {  //upstream?
            quantity<si::velocity, real_t> v_imh = real_t(-.5) * ( // -.5 : averaging + axis orientation
              phc::vterm<real_t>(
                rhod_rr(i,j-1,k) * si::kilograms / si::cubic_metres,
                rhod(i,j-1,k) * si::kilograms / si::cubic_metres,
                rhod(i,rhod.lbound(mtx::j),k) * si::kilograms / si::cubic_metres
              ) +
              phc::vterm<real_t>(
                rhod_rr(i,j,k) * si::kilograms / si::cubic_metres,
                rhod(i,j,k) * si::kilograms / si::cubic_metres,
                rhod(i,rhod.lbound(mtx::j),k) * si::kilograms / si::cubic_metres
              ));
#define upstream_F(a,b,c) ((b)*(c)) // assumes c is negative!
            Rt(i,j,k) -= real_t(
              upstream_F(rhod_rr(i,j,  k), rhod_rr(i,j+1,k), v_iph / si::metres_per_second) -
              upstream_F(rhod_rr(i,j-1,k), rhod_rr(i,j,  k), v_imh / si::metres_per_second)
            ) / real_t(grid.dy() / si::metres);
            v_iph = v_imh;
          }
          Rt(i, rhod.lbound(mtx::j), k) -= real_t(
            rhod_rr(i, rhod.lbound(mtx::j)+1,k)-rhod_rr(i, rhod.lbound(mtx::j),k)
            ) * v_iph / si::metres_per_second / real_t(grid.dy() / si::metres)  ;
// TODO: record accumulated rainfall
        }
      }
    }
  };
*/


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
  : eqs_todo<real_t>(grid, &par), opts(opts), grid(grid)
{ 
  // in activation term all activated droples are assumed to have the radius of 1 um
  const real_t ccnmass = real_t(4./3 * phc::pi<real_t>() * 1e-18) * phc::rho_w<real_t>() / si::kilograms * si::cubic_metres;
  //in autoconversion, ll drizzle drops are assumed to have the radius of 25um
  const real_t drizzle_mass = real_t(4./3 * phc::pi<real_t>() * pow(25.*1e-6, 3)) * phc::rho_w<real_t>() / si::kilograms * si::cubic_metres;

  if (opts[act])  this->sys.at(par.idx_rhod_rv).rhs_terms.push_back(new activation<true>(par, -ccnmass, mean_rd, sdev_rd, n_tot, chem_b));
  if (opts[cond]) this->sys.at(par.idx_rhod_rv).rhs_terms.push_back(new cond_evap<true>(par, -1));
  if (opts[cond] && opts[autoc]) this->sys.at(par.idx_rhod_rv).rhs_terms.push_back(new ce_rain<true>(par, -1));

  if (opts[act])  this->sys.at(par.idx_rhod_th).rhs_terms.push_back(new activation_th(par, ccnmass));
  if (opts[cond]) this->sys.at(par.idx_rhod_th).rhs_terms.push_back(new theta_cond_evap(par));
  if (opts[cond] && opts[autoc]) this->sys.at(par.idx_rhod_th).rhs_terms.push_back(new theta_ce_rain(par));

  this->sys.push_back(new struct eqs<real_t>::gte({
    "rhod_rl", "dry air density times liquid water mixing ratio (i.e. liquid water mass density)",
    this->quan2str(this->par.rho_unit),
    eqs<real_t>::positive_definite(true),
  }));
  if (opts[act])   this->sys.back().rhs_terms.push_back(new activation<false>(par, ccnmass, mean_rd, sdev_rd, n_tot, chem_b));
  if (opts[cond])  this->sys.back().rhs_terms.push_back(new cond_evap<false>(par, 1));
  if (opts[autoc]) this->sys.back().rhs_terms.push_back(new autoconv<true>(par, real_t(-1)));
  if (opts[acc])   this->sys.back().rhs_terms.push_back(new accretion<true>(par, real_t(-1)));
  par.idx_rhod_rl = this->sys.size() - 1;

  this->sys.push_back(new struct eqs<real_t>::gte({
    "rhod_nl", "dry air density times liquid water number concentration",
    this->quan2str(this->par.rho_unit),
    eqs<real_t>::positive_definite(true),
  }));
  if (opts[act])   this->sys.back().rhs_terms.push_back(new activation<false>(par, 1, mean_rd, sdev_rd, n_tot, chem_b));
  if (opts[cond])  this->sys.back().rhs_terms.push_back(new n_evap(par));
  if (opts[autoc] && opts[acc])  this->sys.back().rhs_terms.push_back(new col_sink(par, grid));
  //sink of the clud droplet concentration due to autoconversion is combined with the sink due to accretion
  //TODO have them separately?
  par.idx_rhod_nl = this->sys.size() - 1;                                                                  

  this->sys.push_back(new struct eqs<real_t>::gte({
    "rhod_rr", "dry air density times rain water mixing ratio (i.e. rain water mass density)",
    this->quan2str(this->par.rho_unit),
    eqs<real_t>::positive_definite(true),
  }));
  if (opts[cond] && opts[autoc]) this->sys.back().rhs_terms.push_back(new ce_rain<false>(par, real_t(1)));
  if (opts[autoc]) this-> sys.back().rhs_terms.push_back(new autoconv<false>(par, real_t(1)));
  if (opts[acc])   this-> sys.back().rhs_terms.push_back(new accretion<false>(par, real_t(1)));
//  if (opts[sedi])  this-> sys.back().rhs_terms.push_back(new terminal_velocity(par, grid));
  par.idx_rhod_rr = this->sys.size() - 1;

  this->sys.push_back(new struct eqs<real_t>::gte({
    "rhod_nr", "dry air density times rain water number concentration",
    this->quan2str(this->par.rho_unit),
    eqs<real_t>::positive_definite(true),
  }));
  if (opts[cond]) this->sys.back().rhs_terms.push_back(new nr_evap(par));
  if (opts[autoc]) this-> sys.back().rhs_terms.push_back(new autoconv<false>(par, drizzle_mass));
  par.idx_rhod_nr = this->sys.size() - 1;                                                     

  // auxliary variable for temporary data
  ptr_map_insert(this->aux)("tmp_act", typename eqs<real_t>::axv({
    "tmp_act", "temporary data", "n/a",
    typename eqs<real_t>::invariable(false),
    vector<int>({0, 0, 0})
  })); // TODO: flag it as not to be output

  ptr_map_insert(this->aux)("tmp_cond", typename eqs<real_t>::axv({
    "tmp_cond", "temporary data", "n/a",
    typename eqs<real_t>::invariable(false),
    vector<int>({0, 0, 0})
  })); // TODO: flag it as not to be output

  ptr_map_insert(this->aux)("tmp_ce_rain", typename eqs<real_t>::axv({
    "tmp_ce_rain", "temporary data", "n/a",
    typename eqs<real_t>::invariable(false),
    vector<int>({0, 0, 0})
  })); // TODO: flag it as not to be output

  ptr_map_insert(this->aux)("tmp_auto", typename eqs<real_t>::axv({
    "tmp_auto", "temporary data", "n/a",
    typename eqs<real_t>::invariable(false),
    vector<int>({0, 0, 0})
  })); // TODO: flag it as not to be output

  ptr_map_insert(this->aux)("tmp_acc", typename eqs<real_t>::axv({
    "tmp_acc", "temporary data", "n/a",
    typename eqs<real_t>::invariable(false),
    vector<int>({0, 0, 0})
  })); // TODO: flag it as not to be output

}

// explicit instantiations
#define ICICLE_INSTANTIATE_CLASS eqs_todo_mm
#include "../cmn/cmn_instant.hpp"