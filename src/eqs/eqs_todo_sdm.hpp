/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date April-June 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @brief definition of the eqs_todo_sdm class
 */
#pragma once

#include "../cfg/cfg_boost_odeint.hpp"
#include "../cfg/cfg_thrust.hpp"
#include "eqs_todo.hpp" 
#include "../vel.hpp"
#include "../phc/phc.hpp"
#include "../sdm/sdm_enums.hpp"
#include "../cmn/cmn_error.hpp"

#include <memory>
using std::unique_ptr;

// TODO: option to set the number of threads to use!

/// @brief 
/// implementation of the Super-Droplet Method (@copydetails Shima_et_al_2009, QJRMS, 135)
/// with kappa-Koehler parameterisation of aerosol solubility (@copydetails Petters_and_Kreidenweis_2007, ACP, 7)
/// and ...

template <typename real_t>
class eqs_todo_sdm : public eqs_todo<real_t> 
{
  // a container for storing options (i.e. which processes ar on/off)
  private: typename eqs_todo<real_t>::params par;

  // ctor of eqs_todo_sdm
  public: eqs_todo_sdm(
    const grd<real_t> &grid, 
    const vel<real_t> &velocity,
    map<enum sdm::processes, bool> opts,
    enum sdm::xi_dfntns xi_dfntn,
    enum sdm::ode_algos adve_algo,
    enum sdm::ode_algos sedi_algo,
    enum sdm::ode_algos cond_algo,
    enum sdm::ode_algos chem_algo,
    const int adve_sstp,
    const int sedi_sstp,
    const int cond_sstp,
    const int chem_sstp,
    const int coal_sstp,
    real_t sd_conc_mean,
    real_t min_rd,
    real_t max_rd, 
    real_t mean_rd1, // dry aerosol initial distribution parameters
    real_t mean_rd2, // TODO: encapsulate in a map/enum
    real_t sdev_rd1,
    real_t sdev_rd2,
    real_t n1_tot,
    real_t n2_tot,
    real_t kappa,
    //initial chemical conditions (in air and droplets)
    map<enum sdm::chem_gas, quantity<phc::mixing_ratio, real_t>> opt_gas, 
    map<enum sdm::chem_aq, quantity<si::mass, real_t>> opt_aq 
  )
#if !defined(USE_BOOST_ODEINT) || !defined(USE_THRUST)
  : eqs_todo<real_t>(grid, &this->par)
  {
    error_macro("eqs_todo_sdm requires icicle to be compiled with support for Boost.odeint and Thrust");
  }
#else
  ; 
  // pimpl
  private: struct detail;
  private: unique_ptr<detail> pimpl;
#endif

  void adjustments(
    int n, // TODO: mo¿e jednak bez n...
    vector<ptr_vector<mtx::arr<real_t>>> &psi,
    ptr_unordered_map<string, mtx::arr<real_t>> &aux, 
    const ptr_vector<mtx::arr<real_t>> C,
    const quantity<si::time, real_t> dt,
    bool record
  )
#if !defined(USE_BOOST_ODEINT) || !defined(USE_THRUST)
  { assert(false); }
#else
  ;
#endif
};
