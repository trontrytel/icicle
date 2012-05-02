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

#  include "eqs_todo.hpp"
#  if defined(USE_BOOST_ODEINT)
#    include <boost/numeric/odeint.hpp>
#    include <boost/numeric/odeint/external/thrust/thrust_algebra.hpp>
#    include <boost/numeric/odeint/external/thrust/thrust_operations.hpp>
#    include <boost/numeric/odeint/external/thrust/thrust_resize.hpp>
namespace odeint = boost::numeric::odeint;
#  endif

#  include <thrust/random.h>

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
  typedef thrust::device_vector< value_type > state_type; // TODO: cpu/gpu option

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

  // private fields for ODE
  private: stepper S; 
  private: rhs F;
  private: state_type SD_v; // holds SD parameters that are variable from the ODE point of view (x,y,r,...)

  typedef thrust::device_vector<int>::size_type size_type;

  // private fields for storing SD parameters that are constant from the ODE point of view (n, i, j)
  private: thrust::device_vector<size_type> SD_id;
  private: thrust::device_vector<int> SD_i, SD_j, SD_ij, SD_n;

  // TODO: take from ctor arguments -> command-line-options
  private: const int sd_conc = 1; 
  // TODO: random seed
  

  private: size_type offset_x, offset_y, offset_r;

  // nested functor: RNG
  private: class rng {
    private: thrust::random::taus88 engine; // TODO: RNG engine as an option
    private: thrust::uniform_real_distribution<value_type> dist;
    public: rng(value_type a, value_type b) : dist(a, b) {}
    public: real_t operator()() { return dist(engine); }
  };

  // nested functor: divide by real constant and cast to int
  private: class divide_by_constant
  {
    private: value_type c;
    public: divide_by_constant(value_type c) : c(c) {}
    public: int operator()(value_type x) { return x/c; }
  };

  // nested functor: 
  private: class ravel_indices
  {
    private: int n;
    public: ravel_indices(int n) : n(n) {}
    public: int operator()(int i, int j) { return i + j * n; }
  };

  private: size_type n_part;
  private: const grd<real_t> &grid;

  private: void sd_init()
  {
    // initialise particle positions
    thrust::generate(
      SD_v.begin() + offset_x, 
      SD_v.begin() + offset_x + n_part, 
      rng(0, grid.nx() * (grid.dx() / si::metres))
    ); // x
    thrust::generate(
      SD_v.begin() + offset_y,
      SD_v.begin() + offset_y + n_part, 
      rng(0, grid.ny() * (grid.dy() / si::metres))
    ); // y
    // initialise particle sizes and numbers
    //thrust::generate(SD_v.begin(), SD_v.end(), rng(0, grid.nx() * (grid.dx() / si::metres))); // r
  }

  private: void sd_sync()
  {
    // initialise house-keeping info
    thrust::transform(
      SD_v.begin() + offset_x, 
      SD_v.begin() + offset_x + n_part,
      SD_i.begin(),
      divide_by_constant(grid.dx() / si::metres)
    );
    thrust::transform(
      SD_v.begin() + offset_y, 
      SD_v.begin() + offset_y + n_part,
      SD_j.begin(),
      divide_by_constant(grid.dy() / si::metres)
    );
    thrust::transform(
      SD_i.begin(),
      SD_i.end(),
      SD_j.begin(),
      SD_ij.begin(),
      ravel_indices(grid.ny())
    );
  }
  
  // ctor
  public: eqs_todo_sdm(const grd<real_t> &grid, map<enum processes, bool> opts) :
    eqs_todo<real_t>(grid, &par), opts(opts), grid(grid)
  {
    n_part = grid.nx() * grid.ny() * sd_conc;

    offset_x = 0 * n_part;
    offset_y = 1 * n_part;
    offset_r = 2 * n_part;

    SD_v.resize(3 * n_part);
    SD_i.resize(n_part);
    SD_j.resize(n_part);
    SD_ij.resize(n_part);
    SD_n.resize(n_part);

    // TODO: assert that we use no paralellisation or allow some parallelism!
    // TODO: random seed as an option

    sd_init();

    // auxliary variable for super-droplet conc
    this->aux.push_back(new struct eqs<real_t>::axv({
      "sd_conc", "super-droplet cencentration", this->quan2str(1./si::cubic_metres),
      vector<int>({1, 1, 1}), // dimspan
      vector<int>({0, 0, 0}),  // halo extent
      typename eqs<real_t>::constant(false)
    }));
    //par.idx_dtheta = this->aux.size() - 1;
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

    // 
    sd_sync();

    // velocities?

    // transfer info on thermodynamics to the RHS...

    // do the ODE part (condensation, evaporation, sedimentation)
    S.do_step(boost::ref(F), SD_v, 0, dt / si::seconds);

    // Monte Carlo part (coalescence)

    // retrieve info on thermodynamics from the RHS...
  }

};
#endif
