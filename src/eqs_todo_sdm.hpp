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
// TODO: units in SD_stat, SD_diag, SD_velo... should be possible as odeint supports Boost.units!
// TODO: check why it fails with g++ -O2!!!

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

// @ brief implementation of the Super-Droplet Method by Shima et al. 2009 (QJRMS, 135)
template <typename real_t>
class eqs_todo_sdm : public eqs_todo<real_t> 
{
  // nested class (well... struct)
  protected: struct params : eqs_todo<real_t>::params
  {
    int idx_sd_conc; // auxiliary variables indices
  };

  // a container for storing options (i.e. which processes ar on/off)
  public: enum processes {cond, sedi, coal};

  typedef double value_type; // TODO: option / check if the device supports it
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

  // nested structure: super-droplet velocity field
  struct SD_velo
  {
    // velocity field copy at the device
    thrust::device_vector<value_type> x, y;
    int x_nx, x_ny, y_nx, y_ny; 

    // ctor
    SD_velo(int nx, int ny) 
    {
      x_nx = nx + 1;
      x_ny = ny;
      y_nx = nx;
      y_ny = ny + 1;
      x.resize(x_nx * x_ny);
      y.resize(y_nx * y_ny);
    }
  };

  // nested structure: super-droplet diagnosed variables
  struct SD_diag
  {
    // spectrum moments and corresponding grid cell ids
    thrust::device_vector<int> M_ij;
    thrust::device_vector<size_type> M_0;
    
    // ctor
    SD_diag(int nx, int ny)
    {
      int n_cell = nx * ny;
      M_ij.resize(n_cell);
      M_0.resize(n_cell);
    }
  };

  // nested structure: super-droplet state info
  struct SD_stat
  {
    // number of particles (i.e. super-droplets)
    size_type n_part;

    // SD parameters that are variable from the ODE point of view (x,y,rw,...)
    state_type xy, rw; 

    // ... and since x and y are hidden in one SD.xy, we declare helper iterators
    state_type::iterator x_begin, x_end, y_begin, y_end;

    // SD parameters that are constant from the ODE point of view (n,rd,i,j)
    thrust::device_vector<value_type> rd; // TODO-AJ
    thrust::device_vector<size_type> id, n; // TODO-AJ (n)
    thrust::device_vector<int> i, j, ij;

    // ctor
    SD_stat(int nx, int ny, real_t sd_conc_mean)
    {
      n_part = size_type(real_t(nx * ny) * sd_conc_mean);
      xy.resize(2 * n_part);
      rw.resize(n_part);
      rd.resize(n_part);
      i.resize(n_part);
      j.resize(n_part);
      ij.resize(n_part);
      n.resize(n_part);
      id.resize(n_part);

      x_begin = xy.begin(); 
      x_end   = x_begin + n_part;
      y_begin = x_end;
      y_end   = y_begin + n_part;
    }
  };

  // nested class: RHS of the ODE to be solved to update super-droplet positions
  private: class rhs_xy
  { 
    // nested functor
    private: 
    template <int di, int dj>
    class interpol
    {
      private: const int n;
      private: const state_type &vel;
      private: const SD_stat &stat;
      public: interpol(
        const int n,
        const state_type &vel,
        const SD_stat &stat
      ) : n(n), vel(vel), stat(stat) {}
      public: value_type operator()(size_type id)
      {
        // TODO weighting by position!
        value_type tmp = .5 * (
          *(vel.begin() + (*(stat.i.begin() + id)     ) + (*(stat.j.begin() + id)     ) * n) + 
          *(vel.begin() + (*(stat.i.begin() + id) + di) + (*(stat.j.begin() + id) + dj) * n)
        );
        return tmp;
      }
    };

    // private fields
    private: const SD_stat &stat;
    private: const SD_velo &velo;
 
    // ctor 
    public: rhs_xy(const SD_stat &stat, const SD_velo &velo) : stat(stat), velo(velo) {}

    // overloaded () operator invoked by odeint
    public: void operator()(
      const state_type &xy, state_type &dxy_dt, const value_type
    )
    {
      // TODO use positions to interpolate velocities! (as an option?)
      {
        thrust::counting_iterator<int> iter(0);
        thrust::transform(
          iter, iter + stat.n_part, 
          dxy_dt.begin(), 
          interpol<1,0>(velo.x_nx, velo.x, stat)
        );
      }
      {
        thrust::counting_iterator<int> iter(0);
        thrust::transform(
          iter, iter + stat.n_part, 
          dxy_dt.begin() + stat.n_part, 
          interpol<0,1>(velo.y_nx, velo.y, stat)
        );
      }
    }
  };

  // nested class: RHS of the ODE to be solved to update super-droplet sizes
  private: class rhs_rw
  { 
    // private fields
    private: const SD_stat &stat;

    // ctor
    public: rhs_rw(const SD_stat &stat) : stat(stat) {}
  
    // overloaded () operator invoked by odeint
    public: void operator()(
      const state_type &rw, state_type &drw_dt, const value_type 
    )
    {
      thrust::counting_iterator<int> iter(0);
      thrust::transform(
        iter, iter + stat.n_part, 
        drw_dt.begin() + stat.n_part, 
        drop_growth_equation(stat)
      );
    }
  };

  // nested functor: RNG
  private: class rng {
    // TODO: random seed as an option
    private: thrust::random::taus88 engine; // TODO: RNG engine as an option
    private: thrust::uniform_real_distribution<value_type> dist;
    public: rng(value_type a, value_type b) : dist(a, b) {}
    public: value_type operator()() { return dist(engine); }
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
  private: class modulo
  {
    private: value_type mod;
    public: modulo(value_type mod) : mod(mod) {}
    public: value_type operator()(value_type a) { return fmod(a + mod, mod); }
  };

  // nested functor: 
  private: class copy_from_device 
  {
    private: int n;
    private: const thrust::device_vector<int> &idx2ij;
    private: const thrust::device_vector<size_type> &from;
    private: mtx::arr<real_t> &to;
    public: copy_from_device(int n, 
      const thrust::device_vector<int> &idx2ij,
      const thrust::device_vector<size_type> &from,
      mtx::arr<real_t> &to
    ) : n(n), idx2ij(idx2ij), from(from), to(to) {}
    public: void operator()(int idx) 
    { 
      int ij = *(idx2ij.begin() + idx);
      to(ij % n, ij / n, 0) = *(from.begin() + idx); 
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
//cerr << "to[" << ij << "] = scl * from[" << ij%n << "," << ij/n << ",0]" << endl;
      *(to.begin() + ij) = scl * from(ij % n, ij / n, 0);
    }
  };

  // initialising particle sizes and positions
  private: void sd_init()
  {
    // initialise particle x coordinates
    thrust::generate(
      stat.x_begin, stat.x_end, 
      rng(0, grid.nx() * (grid.dx() / si::metres))
    ); 
    // initialise particle y coordinates
    thrust::generate(
      stat.y_begin, stat.y_end,
      rng(0, grid.ny() * (grid.dy() / si::metres))
    ); 
    // initialise particle dry sizes
    // TODO-AJ
    // initialise particle numbers
    // TODO-AJ
    thrust::fill(stat.n.begin(), stat.n.end(), 1); // TEMP
  }

  // sorting out which particle belongs to which grid cell
  private: void sd_sync()
  {
    // computing stat.i 
    thrust::transform(
      stat.x_begin, stat.x_end,
      stat.i.begin(),
      divide_by_constant(grid.dx() / si::metres)
    );
    // computing stat.j
    thrust::transform(
      stat.y_begin, stat.y_end,
      stat.j.begin(),
      divide_by_constant(grid.dy() / si::metres)
    );
    // computing stat.ij
    thrust::transform(
      stat.i.begin(), stat.i.end(),
      stat.j.begin(),
      stat.ij.begin(),
      ravel_indices(grid.ny())
    );
    // filling-in stat.id
    thrust::sequence(stat.id.begin(), stat.id.end());
    // sorting stat.ij and stat.id
    thrust::sort_by_key(
      stat.ij.begin(), stat.ij.end(),
      stat.id.begin()
    );
    // calculating M0
    thrust::pair<
      thrust::device_vector<int>::iterator, 
      thrust::device_vector<size_type>::iterator
    > n = thrust::reduce_by_key(
      stat.ij.begin(), stat.ij.end(),
      thrust::permutation_iterator<
        thrust::device_vector<size_type>::iterator,
        thrust::device_vector<size_type>::iterator
      >(stat.n.begin(), stat.id.begin()), 
      diag.M_ij.begin(),
      diag.M_0.begin()
    );
  }

  private: void sd_advection(
    const mtx::arr<real_t> &Cx,
    const mtx::arr<real_t> &Cy,
    const quantity<si::time, real_t> dt
  )
  {
    // getting velocities from the Eulerian model (TODO: if constant_velocity -> only once!)
    {
      thrust::counting_iterator<int> iter(0);
      thrust::for_each(iter, iter + velo.x.size(), copy_to_device(
        grid.nx() + 1, Cx, velo.x, (grid.dx() / si::metres) / (dt / si::seconds)
      ));
    }
    {
      thrust::counting_iterator<int> iter(0);
      thrust::for_each(iter, iter + velo.y.size(), copy_to_device(
        grid.nx(), Cy, velo.y, (grid.dx() / si::metres) / (dt / si::seconds)
      ));
    }

    // performing advection using odeint // TODO sedimentation
    S.do_step(boost::ref(F_xy), stat.xy, 0, dt / si::seconds);

    // in-place transforms
    thrust::transform(stat.x_begin, stat.x_end, stat.x_begin, modulo(grid.nx() * (grid.dx() / si::metres)));
    thrust::transform(stat.y_begin, stat.y_end, stat.y_begin, modulo(grid.ny() * (grid.dy() / si::metres)));
  }
  

  // private fields of eqs_todo_sdm
  private: params par;
  private: map<enum processes, bool> opts;
  private: bool constant_velocity;
  private: const grd<real_t> &grid;

  // private fields for ODE machinery
  private: stepper S; 
  private: rhs_xy F_xy;
  private: rhs_rw F_rw;

  // private fields with super droplet structures
  private: SD_stat stat;
  private: SD_velo velo;
  private: SD_diag diag;

  // ctor of eqs_todo_sdm
  public: eqs_todo_sdm(
    const grd<real_t> &grid, 
    const vel<real_t> &velocity,
    map<enum processes, bool> opts,
    real_t sd_conc_mean
  ) : 
    eqs_todo<real_t>(grid, &par), 
    opts(opts), 
    grid(grid), 
    constant_velocity(velocity.is_constant()),
    stat(grid.nx(), grid.ny(), sd_conc_mean),
    velo(grid.nx(), grid.ny()),
    diag(grid.nx(), grid.ny()),
    F_xy(stat, velo),
    F_rw(stat)
  {
    // TODO: assert that we use no paralellisation or allow some parallelism!
    // TODO: random seed as an option
    // TODO: option: number of cores to use (via Thrust)

    // auxliary variable for super-droplet conc
    this->aux.push_back(new struct eqs<real_t>::axv({
      "sd_conc", "super-droplet cencentration", this->quan2str(1./si::cubic_metres),
      typename eqs<real_t>::constant(false),
      vector<int>({0, 0, 0})
    }));
    par.idx_sd_conc = this->aux.size() - 1;

    // initialising super-droplets
    sd_init();
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

    assert(sd_conc.lbound(mtx::k) == sd_conc.ubound(mtx::k)); // 2D

    sd_sync();
    sd_advection(C[0], C[1], dt);

    // foreach below traverses only the grid cells containing super-droplets
    {
      sd_conc(sd_conc.ijk) = real_t(0); // as some grid cells could be void of particles
      thrust::counting_iterator<int> iter(0);
      thrust::for_each(iter, iter + diag.M_ij.size(), copy_from_device(
        grid.nx(), diag.M_ij, diag.M_0, sd_conc // TODO: this is not valid for parallel runs!
      ));
    }


  }

};
#endif
