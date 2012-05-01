/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date April 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef EQS_TODO_SDM_HPP
#  define EQS_TODO_SDM_HPP

// TODO: check if Thrust available
// TODO: something like: #  if !defined(USE_BOOST_ODEINT) error_macro("eqs_todo_sdm requires icicle to be compiled with Boost.odeint");
// TODO: substepping with timesteps as an option

#  include <random> // C++11

#  include "eqs_todo.hpp"
#  if defined(USE_BOOST_ODEINT)
#    include <boost/numeric/odeint.hpp>
#    include <boost/numeric/odeint/external/thrust/thrust_algebra.hpp>
#    include <boost/numeric/odeint/external/thrust/thrust_operations.hpp>
#    include <boost/numeric/odeint/external/thrust/thrust_resize.hpp>
namespace odeint = boost::numeric::odeint;
#  endif

#  include <thrust/random/uniform_real_distribution.h>
#  include <thrust/random/subtract_with_carry_engine.h>

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

  // RHS of the ODE to be solved
  private: class rhs
  { 
    public: void operator()(
      const state_type &x, state_type &dxdt, value_type t
    )
    {
      //F = ...
    }
  };

  // private fields
  private: stepper S; 
  private: rhs F;
  private: state_type X;

  private: const int n_part = 1024; // TODO: as an option
  private: const int n_prop = 3; // x, y, r // TODO: deduce that z is not there...

  private: class rng {
    private: thrust::random::ranlux48_base engine;
    private: thrust::uniform_real_distribution<real_t> dist;
    public: rng(real_t a, real_t b) : dist(a, b) {}
    public: real_t operator()() { return dist(engine); }
  };

  // ctor
  public: eqs_todo_sdm(const grd<real_t> &grid, map<enum processes, bool> opts)
    : eqs_todo<real_t>(grid, &par), opts(opts), X(n_part * n_prop)
  {
    // TODO: assert that we use no paralellisation or allow some parallelism!
    // TODO: random seed as an option

    // initialise positions
    thrust::generate(X.begin()+0*n_part, X.begin()+1*n_part-1, rng(0, grid.nx() * (grid.dx() / si::metres))); // x
    thrust::generate(X.begin()+1*n_part, X.begin()+2*n_part-1, rng(0, grid.ny() * (grid.dy() / si::metres))); // y
    //thrust::generate(X.begin(), X.end(), rng(0, grid.nx() * (grid.dx() / si::metres))); // r

    // TODO: aux variable for density of super-droplets?
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

    // velocities?

    // transfer info on thermodynamics to the RHS...

    // do the ODE part (condensation, evaporation, sedimentation)
    S.do_step(boost::ref(F), X, 0, dt / si::seconds);

    // Monte Carlo part (coalescence)

    // retrieve info on thermodynamics from the RHS...
  }

};
#endif
