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
#  include <thrust/sequence.h>
#  include <thrust/sort.h>

template <typename real_t>
class eqs_todo_sdm : public eqs_todo<real_t> 
{
  // nested class (well... struct)
  protected: struct params : eqs_todo<real_t>::params
  {
    int idx_sd_conc; // auxiliary variables indices
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

  typedef thrust::device_vector<int>::size_type size_type;

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
      //thrust::for_each(dxdt.begin());
    }
  };

  // private fields for ODE machinery
  private: stepper S; 
  private: rhs F;

  // SD parameters that are variable from the ODE point of view (x,y,r,...)
  private: state_type SD_v; 

  // SD parameters that are constant from the ODE point of view (n, i, j)
  private: thrust::device_vector<size_type> SD_id, SD_n;
  private: thrust::device_vector<int> SD_i, SD_j, SD_ij;

  // spectrum moments
  private: thrust::device_vector<int> M_ij;
  private: thrust::device_vector<size_type> M_0;

  // velocity field copy at the device
  private: thrust::device_vector<value_type> U_x, U_y;

  // variables indicating where in SD_v to look for particulat SD properties
  private: size_type offset_x, offset_y, offset_r;

  // nested functor: RNG
  private: class rng {
    // TODO: random seed as an option
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

  // nested functor: multiply by a real constant
  private: class multiply_by_constant
  {
    private: value_type c;
    public: multiply_by_constant(value_type c) : c(c) {}
    public: value_type operator()(value_type x) { return x*c; }
  };

  // nested functor: 
  private: class ravel_indices
  {
    private: int n;
    public: ravel_indices(int n) : n(n) {}
    public: int operator()(int i, int j) { return i + j * n; }
  };

  // nested functor: 
  private: class copy_from_device 
  {
    private: int n;
    private: const thrust::device_vector<size_type> &from;
    private: mtx::arr<real_t> &to;
    public: copy_from_device(int n, 
      const thrust::device_vector<size_type> &from,
      mtx::arr<real_t> &to
    ) : n(n), from(from), to(to) {}
    public: void operator()(int ij) 
    { 
      to(ij % n, ij / n, 0) = *(from.begin() + ij); 
    }
  };

  // nested_functor:
  private: class copy_to_device
  {
    private: int n;
    private: const mtx::arr<real_t> &from;
    private: thrust::device_vector<value_type> &to;
    private: value_type scl;
    public: copy_to_device(int n,
      const mtx::arr<real_t> &from,
      thrust::device_vector<value_type> &to,
      value_type scl = value_type(1)
    ) : n(n), from(from), to(to), scl(scl) {}
    public: void operator()(int ij)
    {
      *(to.begin() + ij) = scl * from(ij % n, ij / n, 0);
    }
  };

  private: size_type n_part;
  private: const grd<real_t> &grid;

  private: void sd_init()
  {
    // initialise particle x coordinates
    thrust::generate(
      SD_v.begin() + offset_x, 
      SD_v.begin() + offset_x + n_part, 
      rng(0, grid.nx() * (grid.dx() / si::metres))
    ); 
    // initialise particle y coordinates
    thrust::generate(
      SD_v.begin() + offset_y,
      SD_v.begin() + offset_y + n_part, 
      rng(0, grid.ny() * (grid.dy() / si::metres))
    ); 
    // initialise particle sizes
    // TODO!
    // initialise particle numbers
    // TODO!
    thrust::fill(SD_n.begin(), SD_n.end(), 1); // TEMP
  }

  private: void sd_sync()
  {
    // computing SD_i 
    thrust::transform(
      SD_v.begin() + offset_x, 
      SD_v.begin() + offset_x + n_part,
      SD_i.begin(),
      divide_by_constant(grid.dx() / si::metres)
    );
    // computing SD_j
    thrust::transform(
      SD_v.begin() + offset_y, 
      SD_v.begin() + offset_y + n_part,
      SD_j.begin(),
      divide_by_constant(grid.dy() / si::metres)
    );
    // computing SD_ij
    thrust::transform(
      SD_i.begin(), SD_i.end(),
      SD_j.begin(),
      SD_ij.begin(),
      ravel_indices(grid.ny())
    );
    // filling-in SD_id
    thrust::sequence(SD_id.begin(), SD_id.end());
    // sorting SD_ij and SD_id
    thrust::sort_by_key(
      SD_ij.begin(), SD_ij.end(),
      SD_id.begin()
    );
    // calculating M0
    thrust::pair<
      thrust::device_vector<int>::iterator, 
      thrust::device_vector<size_type>::iterator
    > n = thrust::reduce_by_key(
      SD_ij.begin(), SD_ij.end(),
      SD_n.begin(),
      M_ij.begin(),
      M_0.begin()
    );
  }
  
  private: bool constant_velocity;

  // ctor
  // TODO: option: number of cores to use (via Thrust)
  public: eqs_todo_sdm(
    const grd<real_t> &grid, 
    const vel<real_t> &velocity,
    map<enum processes, bool> opts,
    int sd_conc_mean
  ) : 
    eqs_todo<real_t>(grid, &par), 
    opts(opts), 
    grid(grid), 
    constant_velocity(velocity.is_constant())
  {
    int n_cell = grid.nx() * grid.ny();
    n_part = n_cell * sd_conc_mean;

    offset_x = 0 * n_part;
    offset_y = 1 * n_part;
    offset_r = 2 * n_part;

    SD_v.resize(3 * n_part);
    SD_i.resize(n_part);
    SD_j.resize(n_part);
    SD_ij.resize(n_part);
    SD_n.resize(n_part);
    SD_id.resize(n_part);

    M_ij.resize(n_cell);
    M_0.resize(n_cell);

    U_x.resize(n_cell + grid.ny());
    U_y.resize(n_cell + grid.nx());

    // TODO: assert that we use no paralellisation or allow some parallelism!
    // TODO: random seed as an option

    sd_init();

    // auxliary variable for super-droplet conc
    this->aux.push_back(new struct eqs<real_t>::axv({
      "sd_conc", "super-droplet cencentration", this->quan2str(1./si::cubic_metres),
      typename eqs<real_t>::constant(false),
      vector<int>({0, 0, 0})
    }));
    par.idx_sd_conc = this->aux.size() - 1;
  }

  public: void adjustments(
    int n, // TODO: mo¿e jednak bez n...
    vector<ptr_vector<mtx::arr<real_t>>> &psi,
    ptr_vector<mtx::arr<real_t>> &aux, 
    const ptr_vector<mtx::arr<real_t>> C,
    const quantity<si::time, real_t> dt
  ) 
  {
    const mtx::arr<real_t>
      &rhod = aux[this->par.idx_rhod];
    mtx::arr<real_t>
      &rhod_th = psi[this->par.idx_rhod_th][n],
      &rhod_rv = psi[this->par.idx_rhod_rv][n],
      //&rhod_rl = psi[this->par.idx_rhod_rl][n],
      //&rhod_rr = psi[this->par.idx_rhod_rr][n],
      &sd_conc = aux[this->par.idx_sd_conc];

    // 
    sd_sync();

    assert(sd_conc.lbound(mtx::k) == sd_conc.ubound(mtx::k)); // 2D

    // foreach below traverses only the grid cells containing super-droplets
    sd_conc(sd_conc.ijk) = real_t(0); 
    thrust::for_each(M_ij.begin(), M_ij.end(), copy_from_device(
      grid.nx(), M_0, sd_conc // TODO: this is not valid for parallel runs!
    ));

    // velocities (TODO: if constant_velocity -> only once!)
    thrust::sequence(U_x.begin(), U_x.end()); // filling with (ij + j)
    thrust::for_each(U_x.begin(), U_x.end(), copy_to_device(
      grid.nx() + 1, C[0], U_x, (dt / si::seconds) / (grid.dx() / si::metres)
    ));
    thrust::sequence(U_y.begin(), U_y.end()); // filling with (ij + j)
    thrust::for_each(U_y.begin(), U_y.end(), copy_to_device(
      grid.nx(), C[1], U_y, (dt / si::seconds) / (grid.dx() / si::metres)
    ));

    // transfer info on thermodynamics to the RHS...

    // do the ODE part (condensation, evaporation, sedimentation)
    S.do_step(boost::ref(F), SD_v, 0, dt / si::seconds);

    // Monte Carlo part (coalescence)

    // retrieve info on thermodynamics from the RHS...
  }

};
#endif
