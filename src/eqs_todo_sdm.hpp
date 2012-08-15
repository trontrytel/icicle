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

#include "eqs_todo.hpp" 
#include "vel.hpp"
#include "phc.hpp"

# if defined(USE_BOOST_ODEINT) && defined(USE_THRUST)
#  include "sdm_functors.hpp"
#  include "sdm_enums.hpp"
#  include "sdm_ode.hpp"
#endif

#include <memory>
using std::unique_ptr;

// TODO: option to set the number of threads to use!

/// @brief 
/// implementation of the Super-Droplet Method (@copydetails Shima_et_al_2009, QJRMS, 135)
/// with kappa-Koehler parameterisation of aerosol solubility (@copydetails Petters_and_Kreidenweis_2007, ACP, 7)
/// and ...

template <typename real_t, int thrust_device_system>
class eqs_todo_sdm : public eqs_todo<real_t> 
{
  // a container for storing options (i.e. which processes ar on/off)
  private: typename eqs_todo<real_t>::params par;
  private: struct detail;

  // ctor of eqs_todo_sdm
  public: eqs_todo_sdm(
    const grd<real_t> &grid, 
    const vel<real_t> &velocity,
    map<enum sdm::processes, bool> opts,
    enum sdm::xi_dfntns xi_dfntn,
    enum sdm::ode_algos xy_algo,
    enum sdm::ode_algos ys_algo,
    enum sdm::ode_algos xi_algo,
    enum sdm::ode_algos chem_algo,
    real_t sd_conc_mean,
    real_t min_rd,
    real_t max_rd, 
    real_t mean_rd1, // dry aerosol initial distribution parameters
    real_t mean_rd2,
    real_t sdev_rd1,
    real_t sdev_rd2,
    real_t n1_tot,
    real_t n2_tot,
    real_t kappa,
    //initial chemical conditions (in air and droplets)
    map<enum sdm::chem_gas, quantity<phc::mixing_ratio, real_t>> opt_gas, 
    map<enum sdm::chem_aq, quantity<divide_typeof_helper<si::amount, si::volume>::type, real_t>> opt_aq 
  )
#if !defined(USE_BOOST_ODEINT) || !defined(USE_THRUST)
: eqs_todo<real_t>(grid, &this->par)
{
  error_macro("eqs_todo_sdm requires icicle to be compiled with support for Boost.odeint and Thrust");
}
#else
; 
#endif

#if defined(USE_BOOST_ODEINT) && defined(USE_THRUST)
  typedef thrust::device_vector<int>::size_type thrust_size_t;

  // private fields of eqs_todo_sdm
  private: map<enum sdm::processes, bool> opts;
  private: bool constant_velocity;
  private: const grd<real_t> &grid;

  // private fields for ODE machinery
  private: unique_ptr<sdm::ode<real_t>> 
    F_adve, F_cond, F_sedi, F_chem;

  // private fields with super droplet structures
  private: sdm::stat_t<real_t> stat;
  private: sdm::envi_t<real_t> envi;

  private: enum sdm::xi_dfntns xi_dfntn;
  private: real_t seed = 1234.;// TODO: option!

  // private field with temporary space
  thrust::device_vector<int> tmp_shrt; // e.g. for grid cell indices
  thrust::device_vector<thrust_size_t> tmp_long;
  thrust::device_vector<real_t> tmp_real;

  public: void adjustments(
    int n, // TODO: mo¿e jednak bez n...
    vector<ptr_vector<mtx::arr<real_t>>> &psi,
    ptr_unordered_map<string, mtx::arr<real_t>> &aux, 
    const ptr_vector<mtx::arr<real_t>> C,
    const quantity<si::time, real_t> dt
  );
#endif
};
