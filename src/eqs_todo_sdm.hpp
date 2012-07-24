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

# if defined(USE_BOOST_ODEINT) && defined(USE_THRUST)
#  include "sdm_functors.hpp"
#  include "sdm_ode.hpp"
#endif

// TODO: option to set the number of threads to use!

/// @brief 
/// implementation of the Super-Droplet Method (@copydetails Shima_et_al_2009, QJRMS, 135)
/// with kappa-Koehler parameterisation of aerosol solubility (@copydetails Petters_and_Kreidenweis_2007, ACP, 7)
/// and ...

namespace sdm 
{
  enum chem_gas {gSO2, gO3, gH2O2};
  enum chem_aq {H, OH, SO2, O3, H2O2, HSO3, SO3};
};

template <typename real_t>
class eqs_todo_sdm : public eqs_todo<real_t> 
{
  // a container for storing options (i.e. which processes ar on/off)
  public: enum processes {adve, cond, sedi, coal, chem};
  public: enum ode_algos {euler, mmid, rk4};
  public: enum xi_dfntns {id, ln, p2, p3};
  private: typename eqs_todo<real_t>::params par;
  private: enum xi_dfntns xi_dfntn;
  private: real_t seed = 1234.;
  private: struct detail;

  // ctor of eqs_todo_sdm
  public: eqs_todo_sdm(
    const grd<real_t> &grid, 
    const vel<real_t> &velocity,
    map<enum processes, bool> opts,
    enum xi_dfntns xi_dfntn,
    enum ode_algos xy_algo,
    enum ode_algos ys_algo,
    enum ode_algos xi_algo,
    enum ode_algos chem_algo,
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

  // initialising particle positions, numbers and dry radii
  private: void sd_init(
    const real_t min_rd, const real_t max_rd,
    const real_t mean_rd1, const real_t mean_rd2,
    const real_t sdev_rd1, const real_t sdev_rd2,
    const real_t n1_tot, const real_t n2_tot, 
    const real_t sd_conc_mean,
    const real_t kappa,
    map<enum sdm::chem_gas, quantity<phc::mixing_ratio, real_t>> opt_gas, 
    map<enum sdm::chem_aq, quantity<divide_typeof_helper<si::amount, si::volume>::type, real_t>> opt_aq
  );

  private: void sd_sync_in(
    const mtx::arr<real_t> &rhod,
    const mtx::arr<real_t> &rhod_th,
    const mtx::arr<real_t> &rhod_rv
  );
  private: void sd_sync_out(
    mtx::arr<real_t> &rhod_th,
    mtx::arr<real_t> &rhod_rv
  );

  // computing diagnostics
  private: void sd_diag(ptr_unordered_map<string, mtx::arr<real_t>> &aux);

  private: void sd_advection(
    const quantity<si::time, real_t> dt,
    const mtx::arr<real_t> &Cx,
    const mtx::arr<real_t> &Cy
  );

  private: void sd_sedimentation(
    const quantity<si::time, real_t> dt,
    const mtx::arr<real_t> &rhod
  );
  
  private: void sd_condevap(
    const quantity<si::time, real_t> dt
  );

  private: void sd_chem(
    const quantity<si::time, real_t> dt
  );

  private: void sd_coalescence(
    const quantity<si::time, real_t> dt
  );

  private: void sd_periodic_x();
  private: void sd_recycle();
  private: void sd_sort();
  private: void sd_shuffle_and_sort();

  // sorting out which particle belongs to which grid cell
  private: void sd_ij();

  // private fields of eqs_todo_sdm
  private: map<enum processes, bool> opts;
  private: bool constant_velocity;
  private: const grd<real_t> &grid;

  // private fields for ODE machinery
  private: unique_ptr<sdm::ode<real_t>> 
    F_adve, F_cond, F_sedi, F_chem;

  // private fields with super droplet structures
  private: sdm::stat_t<real_t> stat;
  private: sdm::envi_t<real_t> envi;

  // private field with temporary space
  thrust::device_vector<int> tmp_shrt; // e.g. for grid cell indices
  thrust::device_vector<thrust_size_t> tmp_long, tmp_long_id; 
  thrust::device_vector<real_t> tmp_real, tmp_real_id; 

  public: void adjustments(
    int n, // TODO: mo¿e jednak bez n...
    vector<ptr_vector<mtx::arr<real_t>>> &psi,
    ptr_unordered_map<string, mtx::arr<real_t>> &aux, 
    const ptr_vector<mtx::arr<real_t>> C,
    const quantity<si::time, real_t> dt
  )
  {
    // TODO: assert(sd_conc.lbound(mtx::k) == sd_conc.ubound(mtx::k)); // 2D
    // TODO: sd_breakup(dt); 
    // TODO: substepping with different timesteps as an option
    // TODO: which order would be best? (i.e. cond, chem, coal, ...)
    const mtx::arr<real_t> &rhod = aux.at("rhod");
    mtx::arr<real_t> &rhod_th = psi[this->par.idx_rhod_th][n];
    mtx::arr<real_t> &rhod_rv = psi[this->par.idx_rhod_rv][n];

    sd_ij();
    sd_sync_in(rhod, rhod_th, rhod_rv);
    if (opts[cond]) sd_condevap(dt); // does init() at first time step - has to be placed after sync, and before others
    if (opts[adve]) 
    {
      sd_advection(dt, C[0], C[1]); 
    }
    if (opts[sedi]) 
    {
      sd_sedimentation(dt, aux.at("rhod")); // TODO: SD recycling!
    }
    sd_periodic_x();
    sd_recycle();
    sd_ij();
    if (opts[chem]) sd_chem(dt);
    if (opts[coal])
    {
      int n_steps = 5;
      for (int i = 0; i < n_steps; ++i)
      {
        sd_shuffle_and_sort();
        sd_coalescence(dt / real_t(n_steps));
      }
    }
    else sd_sort();
    sd_diag(aux); // TODO: only when recording!
    // transfer new rhod_rv and rhod_th back to Blitz
    sd_sync_out(rhod_th, rhod_rv); // has to be placed after condensation?
  }
#endif
};
