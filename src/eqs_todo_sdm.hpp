/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
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
// TODO: units in SD_stat, SD_diag, SD_envi... should be possible as odeint supports Boost.units!
// TODO: check why random number generation fails with g++ -O2!!!

#  include "eqs_todo.hpp"
#  include "phc_lognormal.hpp"

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

// @brief implementation of the Super-Droplet Method by Shima et al. 2009 (QJRMS, 135)
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
  public: enum ode_algos {euler, rk4};
  public: enum xi_dfntns {id, ln, p2, p3};

  typedef double thrust_real_t; // TODO: option / check if the device supports it

  typedef thrust::device_vector<int>::size_type thrust_size_type;

  // nested structure: super-droplet environment (velocity, temperature and moisture field)
  struct SD_envi
  {
    // velocity field (copied from an Arakawa-C grid)
    thrust::device_vector<thrust_real_t> vx, vy;
    int vx_nx, vx_ny, vy_nx, vy_ny; 

    // temperature and water vapour density fields (copied from an Arakawa-C grid)
    thrust::device_vector<thrust_real_t> thta, rhov;

    // ctor
    SD_envi(int nx, int ny) 
    {
      vx_nx = nx + 1;
      vx_ny = ny;
      vy_nx = nx;
      vy_ny = ny + 1;
      vx.resize(vx_nx * vx_ny);
      vy.resize(vy_nx * vy_ny);

      thta.resize(nx * ny); // TODO + 2 x halo if interpolating 
      rhov.resize(nx * ny); // TODO + 2 x halo if interpolating
    }
  };

  // nested structure: super-droplet diagnosed variables
  struct SD_diag
  {
    // spectrum moments and corresponding grid cell ids
    thrust::device_vector<int> M_ij;
    thrust::device_vector<thrust_size_type> M_0;
    
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
    thrust_size_type n_part;

    // SD parameters that are variable from the ODE point of view (x,y,rw,...)
    thrust::device_vector<thrust_real_t> xy, xi; 

    // ... and since x and y are hidden in one SD.xy, we declare helper iterators
    thrust::device_vector<thrust_real_t>::iterator x_begin, x_end, y_begin, y_end;

    // SD parameters that are constant from the ODE point of view (n,rd3,i,j)
    thrust::device_vector<thrust_real_t> rd3; // rd3 -> dry radius to the third power 
    thrust::device_vector<thrust_size_type> id, n; // n -> number of real droplets in a super-droplet
    thrust::device_vector<int> i, j, ij; // location of super-droplet within the grid

    // SD_stat ctor
    SD_stat(int nx, int ny, real_t sd_conc_mean)
    {
      n_part = thrust_size_type(real_t(nx * ny) * sd_conc_mean);
      xy.resize(2 * n_part);
      xi.resize(n_part);
      rd3.resize(n_part);
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

/*
  On May 9, 2012, at 7:44 PM, Karsten Ahnert wrote:
  > ... unfortunately the Rosenbrock method cannot be used with any other state type than ublas.matrix.
  > ... I think, the best steppers for stiff systems and thrust are the
  > runge_kutta_fehlberg78 or the bulirsch_stoer with a very high order. But
  > should benchmark both steppers and choose the faster one.
*/

  typedef odeint::euler<
    thrust::device_vector<thrust_real_t>, // state type
    thrust_real_t, // value_type
    thrust::device_vector<thrust_real_t>, // deriv type
    thrust_real_t, // time type
    odeint::thrust_algebra, 
    odeint::thrust_operations
  > algo_euler;

  typedef odeint::runge_kutta4<
    thrust::device_vector<thrust_real_t>, // state type
    thrust_real_t, // value_type
    thrust::device_vector<thrust_real_t>, // deriv type
    thrust_real_t, // time type
    odeint::thrust_algebra, 
    odeint::thrust_operations
  > algo_rk4;

  // nested class: 
  private: class ode 
  {
    public: virtual void advance(
      thrust::device_vector<thrust_real_t> &x, 
      const quantity<si::time, real_t> &dt
    ) = 0;
  };

  // nested class: 
  private: 
  template <class algo> 
  class ode_algo : public ode
  {
    private: algo stepper;

    // pure virtual method
    public: virtual void operator()(
      const thrust::device_vector<thrust_real_t> &xy, 
      thrust::device_vector<thrust_real_t> &dxy_dt, 
      const thrust_real_t
    ) = 0;

    public: void advance(
      thrust::device_vector<thrust_real_t> &x, 
      const quantity<si::time, real_t> &dt
    )
    {
      stepper.do_step(boost::ref(*this), x, 0, dt / si::seconds);
    }
  };


  // identity: xi = rw
  private: template <typename T = thrust_real_t> struct xi_id 
  { 
    T xi(const T &rw) { return rw; }
    T rw(const T &xi) { return xi; } 
    T dxidrw(const T &xi) { return 1; }
  };

  // xi = ln(rw / nm)
  private: template <typename T = thrust_real_t> struct xi_ln
  { 
    T xi(const T &rw) { return log(rw / T(1e-9)); }
    T rw(const T &xi) { return T(1e-9) * exp(xi); } 
    T dxidrw(const T &rw) { return T(1) / rw; }
  };

  // xi = pow(rw, 2) 
  private: template <typename T = thrust_real_t> struct xi_p2
  { 
    T xi(const T &rw) { return rw * rw; }
    T rw(const T &xi) { return sqrt(xi); } 
    T dxidrw(const T &rw) { return T(2) * rw; }
  };

  // xi = pow(rw, 3) 
  private: template <typename T = thrust_real_t> struct xi_p3
  { 
    T xi(const T &rw) { return rw * rw * rw; }
    T rw(const T &xi) { return pow(xi, T(1)/T(3)); } 
    T dxidrw(const T &rw) { return T(3) * rw * rw; }
  };


  // nested class: the ODE system to be solved to update super-droplet positions
  private: 
  template <class algo>
  class ode_xy : public ode_algo<algo>
  { 
    // nested functor 
    private: 
    template <int di, int dj>
    class interpol
    {
      // private fields
      private: const int n;
      private: const thrust::device_vector<thrust_real_t> &vel;
      private: const SD_stat &stat;

      // ctor
      public: interpol(
        const int n,
        const thrust::device_vector<thrust_real_t> &vel,
        const SD_stat &stat
      ) : n(n), vel(vel), stat(stat) {}

      //  overloaded () operator invoked by thrust::transform()
      public: thrust_real_t operator()(thrust_size_type id)
      {
        // TODO weighting by position!
        return .5 * (
          vel[(stat.i[id]     ) + (stat.j[id]     ) * n] + 
          vel[(stat.i[id] + di) + (stat.j[id] + dj) * n]
        );
      }
    };

    // private fields
    private: const SD_stat &stat;
    private: const SD_envi &envi;
 
    // ctor 
    public: ode_xy(const SD_stat &stat, const SD_envi &envi) : stat(stat), envi(envi) {}

    // overloaded () operator invoked by odeint
    public: void operator()(
      const thrust::device_vector<thrust_real_t>&, 
      thrust::device_vector<thrust_real_t> &dxy_dt, 
      const thrust_real_t
    )
    {
      // TODO use positions to interpolate velocities! (as an option?)
      {
        thrust::counting_iterator<int> iter(0);
        thrust::transform(
          iter, iter + stat.n_part, 
          dxy_dt.begin(), 
          interpol<1,0>(envi.vx_nx, envi.vx, stat)
        );
      }
      {
        thrust::counting_iterator<int> iter(0);
        thrust::transform(
          iter, iter + stat.n_part, 
          dxy_dt.begin() + stat.n_part, 
          interpol<0,1>(envi.vy_nx, envi.vy, stat)
        );
      }
    }
  };

  // nested class: the ODE system to be solved to update super-droplet sizes
  private: 
  template <class algo, class xi>
  class ode_xi : public ode_algo<algo>
  { 
    // nested functor
    private: 
    class drop_growth_equation : public xi
    {
      // private fields
      private: const SD_stat &stat;

      // ctor
      public: drop_growth_equation(const SD_stat &stat) : stat(stat) {}

      // overloaded () operator invoked by thrust::transform()
      public: thrust_real_t operator()(thrust_size_type id)
      {
        return 0;
        //return this->dxidrw(this->rw(stat.xi[id]));
        //  * Maxwell-Mason;
      }
    };

    // private fields
    private: const SD_stat &stat;

    // ctor
    public: ode_xi(const SD_stat &stat) : stat(stat) {}
  
    // overloaded () operator invoked by odeint
    public: void operator()(
      const thrust::device_vector<thrust_real_t>&, 
      thrust::device_vector<thrust_real_t> &dxi_dt, 
      const thrust_real_t 
    )
    {
      thrust::counting_iterator<int> iter(0);
      thrust::transform(
        iter, iter + stat.n_part, 
        dxi_dt.begin(), 
        drop_growth_equation(stat)
      );
    }
  };

  // nested functor: RNG
  private: class rng {
    // TODO: random seed as an option
    private: thrust::random::taus88 engine; // TODO: RNG engine as an option
    private: thrust::uniform_real_distribution<thrust_real_t> dist;
    public: rng(thrust_real_t a, thrust_real_t b) : dist(a, b) {}
    public: thrust_real_t operator()() { return dist(engine); }
  };

  // nested functor: divide by real constant and cast to int
  private: class divide_by_constant
  {
    private: thrust_real_t c;
    public: divide_by_constant(thrust_real_t c) : c(c) {}
    public: int operator()(thrust_real_t x) { return x/c; }
  };

  // nested functor: multiply by a real constant
  private: class multiply_by_constant
  {
    private: thrust_real_t c;
    public: multiply_by_constant(thrust_real_t c) : c(c) {}
    public: thrust_real_t operator()(thrust_real_t x) { return x*c; }
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
    private: thrust_real_t mod;
    public: modulo(thrust_real_t mod) : mod(mod) {}
    public: thrust_real_t operator()(thrust_real_t a) { return fmod(a + mod, mod); }
  };

  // nested functor: 
  private: class copy_from_device 
  {
    private: int n;
    private: const thrust::device_vector<int> &idx2ij;
    private: const thrust::device_vector<thrust_size_type> &from;
    private: mtx::arr<real_t> &to;
    private: real_t scl;

    // ctor
    public: copy_from_device(int n, 
      const thrust::device_vector<int> &idx2ij,
      const thrust::device_vector<thrust_size_type> &from,
      mtx::arr<real_t> &to,
      real_t scl = real_t(1)
    ) : n(n), idx2ij(idx2ij), from(from), to(to), scl(scl) {}

    public: void operator()(int idx) 
    { 
      to(idx2ij[idx] % n, idx2ij[idx] / n, 0) = scl * from[idx]; 
    }
  };

  // nested_functor:
  private: class copy_to_device
  {
    private: int n;
    private: const mtx::arr<real_t> &from;
    private: thrust::device_vector<thrust_real_t> &to;
    private: thrust_real_t scl;

    // ctor
    public: copy_to_device(int n,
      const mtx::arr<real_t> &from,
      thrust::device_vector<thrust_real_t> &to,
      thrust_real_t scl = thrust_real_t(1)
    ) : n(n), from(from), to(to), scl(scl) {}

    public: void operator()(int ij)
    {
      to[ij] = scl * from(ij % n, ij / n, 0);
    }
  };

  // nested functor
  private: class lognormal
  {
    private: quantity<si::length, real_t> mean_rd1, mean_rd2;
    private: quantity<si::dimensionless, real_t> sdev_rd1, sdev_rd2;
    private: quantity<power_typeof_helper<si::length, static_rational<-3>>::type, real_t> n1_tot, n2_tot;
    
    public: lognormal(
      const quantity<si::length, real_t> mean_rd1,
      const quantity<si::dimensionless, real_t> sdev_rd1,
      const quantity<power_typeof_helper<si::length, static_rational<-3>>::type, real_t> n1_tot,
      const quantity<si::length, real_t> mean_rd2,
      const quantity<si::dimensionless, real_t> sdev_rd2,
      const quantity<power_typeof_helper<si::length, static_rational<-3>>::type, real_t> n2_tot
    ) : mean_rd1(mean_rd1), sdev_rd1(sdev_rd1), n1_tot(n1_tot), mean_rd2(mean_rd2), sdev_rd2(sdev_rd2), n2_tot(n2_tot) {}

    public: thrust_real_t operator()(real_t lnrd)
    {
      return thrust_real_t(( //TODO allow more modes of distribution
// TODO: logarith or not: as an option
          phc::log_norm_n_e<real_t>(mean_rd1, sdev_rd1, n1_tot, lnrd) + 
          phc::log_norm_n_e<real_t>(mean_rd2, sdev_rd2, n2_tot, lnrd) 
        ) * si::cubic_metres
      );
    }
  };

  // initialising particle sizes and positions
  private: void sd_init(
    real_t min_rd, real_t max_rd,
    real_t mean_rd1, real_t mean_rd2,
    real_t sdev_rd1, real_t sdev_rd2,
    real_t n1_tot, real_t n2_tot, 
    real_t sd_conc_mean)
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
    // initialise particle dry sizes (temporarily with logarithms of radius!)
    thrust::generate(
      stat.rd3.begin(), stat.rd3.end(), 
      rng(log(min_rd), log(max_rd))
    ); 
    // initialise particle numbers  
    // TODO: assert that the distribution is < epsilon at rd_min and rd_max
    real_t multi = log(max_rd/min_rd) / sd_conc_mean 
      * (grid.dx() * grid.dy() * grid.dz() / si::cubic_metres); 
    thrust::transform(
      stat.rd3.begin(), stat.rd3.end(), 
      stat.n.begin(), 
      lognormal(
        mean_rd1 * si::metres, sdev_rd1, multi * n1_tot / si::cubic_metres, 
        mean_rd2 * si::metres, sdev_rd2, multi * n2_tot / si::cubic_metres
      ) 
    );
    // converting rd back from logarihms to rd3 
    {
      struct exp3x { thrust_real_t operator()(thrust_real_t x) { return exp(3*x); } };
      thrust::transform(stat.rd3.begin(), stat.rd3.end(), stat.rd3.begin(), exp3x());
    }
    // initialise particle wet radii (the xi variable)
    // TODO!
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
      thrust::device_vector<thrust_size_type>::iterator
    > n = thrust::reduce_by_key(
      stat.ij.begin(), stat.ij.end(),
      thrust::permutation_iterator<
        thrust::device_vector<thrust_size_type>::iterator,
        thrust::device_vector<thrust_size_type>::iterator
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
      thrust::for_each(iter, iter + envi.vx.size(), copy_to_device(
        grid.nx() + 1, Cx, envi.vx, (grid.dx() / si::metres) / (dt / si::seconds)
      ));
    }
    {
      thrust::counting_iterator<int> iter(0);
      thrust::for_each(iter, iter + envi.vy.size(), copy_to_device(
        grid.nx(), Cy, envi.vy, (grid.dx() / si::metres) / (dt / si::seconds)
      ));
    }

    // performing advection using odeint // TODO sedimentation
    F_xy->advance(stat.xy, dt);

    // periodic boundary conditions (in-place transforms)
    thrust::transform(stat.x_begin, stat.x_end, stat.x_begin, modulo(grid.nx() * (grid.dx() / si::metres)));
    thrust::transform(stat.y_begin, stat.y_end, stat.y_begin, modulo(grid.ny() * (grid.dy() / si::metres)));
  }
  
  private: void sd_condevap(
    const mtx::arr<real_t> &thta,
    const mtx::arr<real_t> &rhov,
    const quantity<si::time, real_t> dt
  )
  {
    // getting thermodynamic fields from the Eulerian model 
    {
      thrust::counting_iterator<int> iter(0);
      thrust::for_each(iter, iter + envi.thta.size(), copy_to_device(
        grid.nx(), thta, envi.thta
      ));
    }
    {
      thrust::counting_iterator<int> iter(0);
      thrust::for_each(iter, iter + envi.rhov.size(), copy_to_device(
        grid.nx(), rhov, envi.rhov
      ));
    }

    // growing/shrinking the droplets
    F_xi->advance(stat.xi, dt);
  }

  // private fields of eqs_todo_sdm
  private: params par;
  private: map<enum processes, bool> opts;
  private: bool constant_velocity;
  private: const grd<real_t> &grid;
  private: real_t min_rd;
  private: real_t max_rd;
  private: real_t mean_rd1;  // dry aerosol initial distribution parameters
  private: real_t mean_rd2;
  private: real_t sdev_rd1;
  private: real_t sdev_rd2;
  private: real_t n1_tot;
  private: real_t n2_tot;
  private: real_t sd_conc_mean;

  // private fields for ODE machinery
  private: unique_ptr<ode> F_xy;
  private: unique_ptr<ode> F_xi; 

  // private fields with super droplet structures
  private: SD_stat stat;
  private: SD_envi envi;
  private: SD_diag diag;

  // ctor of eqs_todo_sdm
  public: eqs_todo_sdm(
    const grd<real_t> &grid, 
    const vel<real_t> &velocity,
    map<enum processes, bool> opts,
    enum xi_dfntns xi_dfntn,
    enum ode_algos xy_algo,
    enum ode_algos xi_algo,
    real_t sd_conc_mean,
    real_t min_rd,
    real_t max_rd, 
    real_t mean_rd1, // dry aerosol initial distribution parameters
    real_t mean_rd2,
    real_t sdev_rd1,
    real_t sdev_rd2,
    real_t n1_tot,
    real_t n2_tot
  ) :
    eqs_todo<real_t>(grid, &par), 
    opts(opts), 
    grid(grid), 
    constant_velocity(velocity.is_constant()),
    stat(grid.nx(), grid.ny(), sd_conc_mean),
    envi(grid.nx(), grid.ny()),
    diag(grid.nx(), grid.ny()),
    min_rd(min_rd), max_rd(max_rd),
    mean_rd1(mean_rd1), mean_rd2(mean_rd2),
    sdev_rd1(sdev_rd1), sdev_rd2(sdev_rd2),
    n1_tot(n1_tot), n2_tot(n2_tot),
    sd_conc_mean(sd_conc_mean)
  {
    // TODO: assert that we use no paralellisation or allow some parallelism!
    // TODO: random seed as an option
    // TODO: option: number of cores to use (via Thrust)

    // auxliary variable for super-droplet conc
    this->aux.push_back(new struct eqs<real_t>::axv({
      "sd_conc", "particle cencentration", this->quan2str(real_t(1)/si::cubic_metres),
      typename eqs<real_t>::constant(false),
      vector<int>({0, 0, 0})
    }));
    par.idx_sd_conc = this->aux.size() - 1;

    switch (xy_algo)
    {
      case euler: F_xy.reset(new ode_xy<algo_euler>(stat, envi)); break;
      case rk4  : F_xy.reset(new ode_xy<algo_rk4  >(stat, envi)); break;
      default: assert(false);
    }
    switch (xy_algo)
    {
      case euler: switch (xi_dfntn)
      {
        case id : F_xi.reset(new ode_xi<algo_euler, xi_id<>>(stat)); break;
        case ln : F_xi.reset(new ode_xi<algo_euler, xi_ln<>>(stat)); break;
        case p2 : F_xi.reset(new ode_xi<algo_euler, xi_p2<>>(stat)); break;
        case p3 : F_xi.reset(new ode_xi<algo_euler, xi_p3<>>(stat)); break;
        default: assert(false);
      }
      case rk4  : switch (xi_dfntn) 
      {
        case id : F_xi.reset(new ode_xi<algo_rk4,   xi_id<>>(stat)); break;
        case ln : F_xi.reset(new ode_xi<algo_rk4,   xi_ln<>>(stat)); break;
        case p2 : F_xi.reset(new ode_xi<algo_rk4,   xi_p2<>>(stat)); break;
        case p3 : F_xi.reset(new ode_xi<algo_rk4,   xi_p3<>>(stat)); break;
        default: assert(false);
      }
      default: assert(false);
    } 

    // initialising super-droplets
    sd_init(min_rd, max_rd, mean_rd1, mean_rd2, sdev_rd1, sdev_rd2, n1_tot, n2_tot, sd_conc_mean);

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

// TODO: which order would be best?
    sd_sync();
    sd_advection(C[0], C[1], dt); // TODO: sedimentation
    sd_condevap(rhod_th, rhod_rv, dt);
//    sd_coalescence(dt);
//    sd_breakup(dt);

    // TODO: this needs to be done only during record steps... how to handle it?
    // foreach below traverses only the grid cells containing super-droplets
    {
      sd_conc(sd_conc.ijk) = real_t(0); // as some grid cells could be void of particles
      thrust::counting_iterator<int> iter(0);
      thrust::for_each(iter, iter + diag.M_ij.size(), copy_from_device(
        grid.nx(), diag.M_ij, diag.M_0, sd_conc, 
        real_t(1) / grid.dx() / grid.dy() / grid.dz() * si::cubic_metres // TODO: this is not valid for parallel runs!
      ));
    }


  }

};
#endif
