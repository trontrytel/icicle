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
#  include "sdm_ode_xi.hpp"
#  include "sdm_ode_xy.hpp"
#endif

// TODO: option to set the number of threads to use!

/// @brief 
/// implementation of the Super-Droplet Method (@copydetails Shima_et_al_2009, QJRMS, 135)
/// with kappa-Koehler parameterisation of aerosol solubility (@copydetails Petters_and_Kreidenweis_2007, ACP, 7)
/// and ...
template <typename real_t>
class eqs_todo_sdm : public eqs_todo<real_t> 
{
  // a container for storing options (i.e. which processes ar on/off)
  public: enum processes {adve, cond, sedi, coal};
  public: enum ode_algos {euler, rk4};
  public: enum xi_dfntns {id, ln, p2, p3};

  private: typename eqs_todo<real_t>::params par;

  // ctor of eqs_todo_sdm
  public: eqs_todo_sdm(
    const grd<real_t> &grid, 
    const vel<real_t> &velocity,
    map<enum processes, bool> opts,
    enum xi_dfntns xi_dfntn,
    enum ode_algos xy_algo,
    enum ode_algos ys_algo,
    enum ode_algos xi_algo,
    real_t sd_conc_mean,
    real_t min_rd,
    real_t max_rd, 
    real_t mean_rd1, // dry aerosol initial distribution parameters
    real_t mean_rd2,
    real_t sdev_rd1,
    real_t sdev_rd2,
    real_t n1_tot,
    real_t n2_tot,
    real_t kappa
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
    const real_t kappa
  );

  // sorting out which particle belongs to which grid cell
  private: void sd_sync(
    const mtx::arr<real_t> &rhod,
    const mtx::arr<real_t> &rhod_th,
    const mtx::arr<real_t> &rhod_rv
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

  // private fields of eqs_todo_sdm
  private: map<enum processes, bool> opts;
  private: bool constant_velocity;
  private: const grd<real_t> &grid;

  // private fields for ODE machinery
  private: unique_ptr<sdm::ode<real_t>> F_xy;
  private: unique_ptr<sdm::ode<real_t>> F_xi; 
  private: unique_ptr<sdm::ode<real_t>> F_ys; 

  // private fields with super droplet structures
  private: sdm::stat_t<real_t> stat;
  private: sdm::envi_t<real_t> envi;

  // private field with temporary space
  thrust::device_vector<int> tmp_shrt; // e.g. for grid cell indices
  thrust::device_vector<thrust_size_t> tmp_long; // e.g. for particle concentrations

  public: void adjustments(
    int n, // TODO: mo¿e jednak bez n...
    vector<ptr_vector<mtx::arr<real_t>>> &psi,
    ptr_unordered_map<string, mtx::arr<real_t>> &aux, 
    const ptr_vector<mtx::arr<real_t>> C,
    const quantity<si::time, real_t> dt
  )
  {
    //assert(sd_conc.lbound(mtx::k) == sd_conc.ubound(mtx::k)); // 2D

    // TODO: substepping with different timesteps as an option
    // TODO: which order would be best?
    sd_sync(
      aux.at("rhod"),
      psi[this->par.idx_rhod_th][n],
      psi[this->par.idx_rhod_rv][n]
    );
    if (opts[cond]) sd_condevap(dt); // does init() at first time step - has to be placed after sync, and before others
    if (opts[adve]) sd_advection(dt, C[0], C[1]); // includes periodic boundary for super droplets!
    if (opts[sedi]) sd_sedimentation(dt, aux.at("rhod")); // TODO: SD recycling!
//    sd_coalescence(dt);
//    sd_breakup(dt);
    sd_diag(aux); 
  }
#endif
};
