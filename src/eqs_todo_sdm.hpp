/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date April 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef EQS_TODO_SDM_HPP
#  define EQS_TODO_SDM_HPP

#  include "eqs_todo.hpp"
#  if defined(USE_BOOST_ODEINT)
#    include <boost/numeric/odeint.hpp>
#      include <boost/numeric/odeint/external/thrust/thrust_algebra.hpp>
#      include <boost/numeric/odeint/external/thrust/thrust_operations.hpp>
#      include <boost/numeric/odeint/external/thrust/thrust_resize.hpp>
namespace odeint = boost::numeric::odeint;
#  endif

template <typename real_t>
class eqs_todo_sdm : public eqs_todo<real_t> 
{
  // nested class (well... struct)
  protected: struct params : eqs_todo<real_t>::params
  {
    int idx_rhod_rl, idx_rhod_rr; // auxiliary variables indices
  };
  protected: params par;

  // a container for storing options (i.e. which processes ar on/off)
  public: enum processes {cond, sedi, coal};
  private: map<enum processes, bool> opts;

  // ctor
  public: eqs_todo_sdm(const grd<real_t> &grid, map<enum processes, bool> opts)
    : eqs_todo<real_t>(grid, &par), opts(opts)
  {
/*
    this->aux.push_back(new struct eqs<real_t>::gte({
      "rhod_rl", "dry air density times liquid water mixing ratio (i.e. liquid water mass density)",
      this->quan2str(this->par.rho_unit),
      typename eqs<real_t>::positive_definite(true),
    }));
    par.idx_rhod_rl = this->sys.size() - 1;

    this->aux.push_back(new struct eqs<real_t>::gte({
      "rhod_rr", "dry air density times rain water mixing ratio (i.e. rain water mass density)",
      this->quan2str(par.rho_unit),
      typename eqs<real_t>::positive_definite(true),
    }));
    par.idx_rhod_rr = this->sys.size() - 1;
*/
  }

  // RHS of the ODE to be solved for Super-Droplet Evolution
  private: class rhs
  { 
    public: void operator()(
      const quantity<multiply_typeof_helper<si::mass_density, si::temperature>::type, real_t> rhod_th, 
      quantity<si::temperature, real_t> &F, 
      const quantity<si::mass_density, real_t> rhod_rv
    )
    {
      //F = ...
    }
  };

  typedef double value_type; // TODO: option / check if the device supports it
  //typedef thrust::device_vector< value_type > state_type; 
  typedef thrust::host_vector< value_type > state_type; // TODO: cpu/gpu option

  typedef state_type deriv_type;
  typedef value_type time_type;

  //typedef odeint::euler< // TODO: option
  typedef odeint::runge_kutta4<
    state_type,
    value_type,
    deriv_type,
    time_type,
    odeint::thrust_algebra, 
    odeint::thrust_operations
  > stepper;

  public: void adjustments(
    int n, // TODO: mo¿e jednak bez n...
    const ptr_vector<mtx::arr<real_t>> &aux, 
    vector<ptr_vector<mtx::arr<real_t>>> &psi,
    const quantity<si::time, real_t> dt
  ) 
  {
    const mtx::arr<real_t>
      &rhod = aux[this->par.idx_rhod];
    mtx::arr<real_t>
      &rhod_th = psi[this->par.idx_rhod_th][n],
      &rhod_rv = psi[this->par.idx_rhod_rv][n],
      &rhod_rl = psi[this->par.idx_rhod_rl][n],
      &rhod_rr = psi[this->par.idx_rhod_rr][n];

    // TODO: substepping with timesteps as an option
    ode_part(rhod, rhod_th, rhod_rv, rhod_rl, rhod_rr, dt); // condensation/evaporation and sedimentation
    monte_carlo_part(); // coalescence
  }
 
  private: void ode_part(
    const mtx::arr<real_t> &rhod,
    mtx::arr<real_t> &rhod_th,
    mtx::arr<real_t> &rhod_rv,
    mtx::arr<real_t> &rhod_rl,
    mtx::arr<real_t> &rhod_rr,
    const quantity<si::time, real_t> dt
  )   
  {
#  if !defined(USE_BOOST_ODEINT)
    error_macro("eqs_todo_sdm requires icicle to be compiled with Boost.odeint");
#  else

    stepper S; // TODO: would be better to instantiate in the ctor (but what about thread safety! :()
    rhs F;

/*
    S.do_step(
      boost::ref(F), 
      tmp,
      rhod_rv(i,j,k) * si::kilograms / si::cubic_metres, 
      drho_rv        * si::kilograms / si::cubic_metres
    );
*/

#  endif
  }

  private: void monte_carlo_part()
  {
    // TODO...
  }
};
#endif
