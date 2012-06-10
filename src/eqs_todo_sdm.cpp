/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date April-June 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @brief definition of the eqs_todo_sdm class
 */

#include "cfg.hpp"

#if defined(USE_BOOST_ODEINT) && defined(USE_THRUST)

#  include "eqs_todo_sdm.hpp"
#  include "phc_lognormal.hpp" // TODO: not here?
#  include "phc_kappa_koehler.hpp" // TODO: not here?

template <typename real_t>
eqs_todo_sdm<real_t>::eqs_todo_sdm(
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
  real_t n2_tot,
  real_t kappa
) : eqs_todo<real_t>(grid, &this->par), 
  opts(opts), 
  grid(grid), 
  constant_velocity(velocity.is_constant()),
  stat(grid.nx(), grid.ny(), sd_conc_mean),
  envi(grid.nx(), grid.ny())
{
  // TODO: assert that we use no paralellisation or allow some parallelism!
  // TODO: random seed as an option
  // TODO: option: number of cores to use (via Thrust)

  tmp_shrt.resize(grid.nx() * grid.ny());
  tmp_long.resize(grid.nx() * grid.ny());
  tmp_real.resize(grid.nx() * grid.ny());

  // auxliary variable for super-droplet conc
  ptr_map_insert(this->aux)("sd_conc", typename eqs<real_t>::axv({
    "sd_conc", "super-droplet cencentration", this->quan2str(real_t(1)/grid.dx()/grid.dy()/grid.dz()),
    typename eqs<real_t>::invariable(false),
    vector<int>({0, 0, 0})
  }));

  // auxliary variable for total particle concentration
  ptr_map_insert(this->aux)("n_tot", typename eqs<real_t>::axv({
    "n_tot", "total particle cencentration", this->quan2str(real_t(1)/si::cubic_metres),
    typename eqs<real_t>::invariable(false),
    vector<int>({0, 0, 0})
  }));

  // auxliary variable for CCN concentration (temporary kludge)
  ptr_map_insert(this->aux)("n_ccn", typename eqs<real_t>::axv({
    "n_ccn", "CCN cencentration", this->quan2str(real_t(1)/si::cubic_metres),
    typename eqs<real_t>::invariable(false),
    vector<int>({0, 0, 0})
  }));

  // initialising super-droplets
  sd_init(min_rd, max_rd, mean_rd1, mean_rd2, sdev_rd1, sdev_rd2, n1_tot, n2_tot, sd_conc_mean, kappa);

  // initialising ODE right-hand-sides
  switch (xy_algo)
  {
    case euler: F_xy.reset(new sdm::ode_xy<real_t, sdm::algo_euler>(stat, envi)); break;
    case rk4  : F_xy.reset(new sdm::ode_xy<real_t, sdm::algo_rk4  >(stat, envi)); break;
    default: assert(false);
  }
  switch (xy_algo)
  {
    case euler: switch (xi_dfntn)
    {
      case id : F_xi.reset(new sdm::ode_xi<real_t, sdm::algo_euler, sdm::xi_id<>>(stat, envi, tmp_real)); break;
      case ln : F_xi.reset(new sdm::ode_xi<real_t, sdm::algo_euler, sdm::xi_ln<>>(stat, envi, tmp_real)); break;
      case p2 : F_xi.reset(new sdm::ode_xi<real_t, sdm::algo_euler, sdm::xi_p2<>>(stat, envi, tmp_real)); break;
      case p3 : F_xi.reset(new sdm::ode_xi<real_t, sdm::algo_euler, sdm::xi_p3<>>(stat, envi, tmp_real)); break;
      default: assert(false);
    } 
    break;
    case rk4  : switch (xi_dfntn) 
    {
      case id : F_xi.reset(new sdm::ode_xi<real_t, sdm::algo_rk4,   sdm::xi_id<>>(stat, envi, tmp_real)); break;
      case ln : F_xi.reset(new sdm::ode_xi<real_t, sdm::algo_rk4,   sdm::xi_ln<>>(stat, envi, tmp_real)); break;
      case p2 : F_xi.reset(new sdm::ode_xi<real_t, sdm::algo_rk4,   sdm::xi_p2<>>(stat, envi, tmp_real)); break;
      case p3 : F_xi.reset(new sdm::ode_xi<real_t, sdm::algo_rk4,   sdm::xi_p3<>>(stat, envi, tmp_real)); break;
      default: assert(false);
    }
    break;
    default: assert(false);
  } 
}

// initialising particle positions, numbers and dry radii
template <typename real_t>
void eqs_todo_sdm<real_t>::sd_init(
  const real_t min_rd, const real_t max_rd,
  const real_t mean_rd1, const real_t mean_rd2,
  const real_t sdev_rd1, const real_t sdev_rd2,
  const real_t n1_tot, const real_t n2_tot, 
  const real_t sd_conc_mean,
  const real_t kappa
)
{
  thrust_real_t seed = 1234.;

  // initialise particle coordinates
  {
    sdm::rng<thrust_real_t> rand_x(0, grid.nx() * (grid.dx() / si::metres), seed);
    thrust::generate(stat.x_begin, stat.x_end, rand_x);
    seed = rand_x();
  }
  {
    sdm::rng<thrust_real_t> rand_y(0, grid.ny() * (grid.dy() / si::metres), seed);
    thrust::generate(stat.y_begin, stat.y_end, rand_y);
    seed = rand_y();
  }

  // initialise particle dry size spectrum 
  // TODO: assert that the distribution is < epsilon at rd_min and rd_max
  thrust::generate( // rd3 temporarily means logarithms of radius!
    stat.rd3.begin(), stat.rd3.end(), 
    sdm::rng<thrust_real_t>(log(min_rd), log(max_rd), seed) 
  ); 
  {
    real_t multi = log(max_rd/min_rd) / sd_conc_mean 
      * (grid.dx() * grid.dy() * grid.dz() / si::cubic_metres); 
    thrust::transform(
      stat.rd3.begin(), stat.rd3.end(), 
      stat.n.begin(), 
      sdm::lognormal<thrust_real_t>(
        mean_rd1 * si::metres, sdev_rd1, multi * n1_tot / si::cubic_metres, 
        mean_rd2 * si::metres, sdev_rd2, multi * n2_tot / si::cubic_metres
      ) 
    );
    // converting rd back from logarihms to rd3 
    struct exp3x { thrust_real_t operator()(thrust_real_t x) { return exp(3*x); } };
    thrust::transform(stat.rd3.begin(), stat.rd3.end(), stat.rd3.begin(), exp3x());
  }

  // initialise kappas
  thrust::fill(stat.kpa.begin(), stat.kpa.end(), kappa);
}

  // sorting out which particle belongs to which grid cell
template <typename real_t>
void eqs_todo_sdm<real_t>::sd_sync(
  const mtx::arr<real_t> &rhod,
  const mtx::arr<real_t> &rhod_th,
  const mtx::arr<real_t> &rhod_rv
)
{
  // sorting the particles
  // computing stat.i 
  thrust::transform(
    stat.x_begin, stat.x_end,
    stat.i.begin(),
    sdm::divide_by_constant<thrust_real_t>(grid.dx() / si::metres)
  );
  // computing stat.j
  thrust::transform(
    stat.y_begin, stat.y_end,
    stat.j.begin(),
    sdm::divide_by_constant<thrust_real_t>(grid.dy() / si::metres)
  );
  // computing stat.ij
  thrust::transform(
    stat.i.begin(), stat.i.end(),
    stat.j.begin(),
    stat.ij.begin(),
    sdm::ravel_indices(grid.ny())
  );
  // filling-in stat.id
  thrust::sequence(stat.id.begin(), stat.id.end());
  // sorting stat.ij and stat.id
  thrust::sort_by_key(
    stat.ij.begin(), stat.ij.end(),
    stat.id.begin()
  );

  // getting thermodynamic fields from the Eulerian model 
  {
    thrust::counting_iterator<thrust_size_t> iter(0);
    thrust::for_each(iter, iter + envi.rhod.size(), 
      sdm::copy_to_device<real_t, thrust_real_t>(
        grid.nx(), rhod, envi.rhod
      )
    ); // TODO: this could be done just once in the kinematic model!
  }
  {
    thrust::counting_iterator<thrust_size_t> iter(0);
    thrust::for_each(iter, iter + envi.rhod_th.size(), 
      sdm::copy_to_device<real_t, thrust_real_t>(
        grid.nx(), rhod_th, envi.rhod_th
      )
    );
  }
  {
    thrust::counting_iterator<thrust_size_t> iter(0);
    thrust::for_each(iter, iter + envi.rhod_rv.size(), 
      sdm::copy_to_device<real_t, thrust_real_t>(
        grid.nx(), rhod_rv, envi.rhod_rv
      )
    );
  }
}

// computing diagnostics
template <typename real_t>
void eqs_todo_sdm<real_t>::sd_diag(ptr_unordered_map<string, mtx::arr<real_t>> &aux)
{
  // calculating super-droplet concentration (per grid cell)
  {
    // zeroing the temporary var 
    thrust::fill(tmp_long.begin(), tmp_long.end(), 0); // TODO: is it needed?
    // doing the reduction
    thrust::pair<
      thrust::device_vector<int>::iterator, 
      thrust::device_vector<thrust_size_t>::iterator
    > n = thrust::reduce_by_key(
      stat.ij.begin(), stat.ij.end(),
      thrust::make_constant_iterator(1),
      tmp_shrt.begin(), // will store the grid cell indices
      tmp_long.begin()  // will store the concentrations per grid cell
    );
    // writing to aux
    mtx::arr<real_t> &sd_conc = aux.at("sd_conc");
    sd_conc(sd_conc.ijk) = real_t(0); // as some grid cells could be void of particles
    thrust::counting_iterator<thrust_size_t> iter(0);
    thrust::for_each(iter, iter + (n.first - tmp_shrt.begin()), 
      sdm::copy_from_device<real_t>(
        grid.nx(), tmp_shrt, tmp_long, sd_conc 
      )
    );
  }

  // calculating the zero-th moment (i.e. total particle concentration per unit volume)
  {
    // zeroing the temporary var
    thrust::fill(tmp_long.begin(), tmp_long.end(), 0); // TODO: is it needed?
    // doing the reduction
    thrust::pair<
      thrust::device_vector<int>::iterator, 
      thrust::device_vector<thrust_size_t>::iterator
    > n = thrust::reduce_by_key(
      stat.ij.begin(), stat.ij.end(),
      thrust::permutation_iterator<
        thrust::device_vector<thrust_size_t>::iterator,
        thrust::device_vector<thrust_size_t>::iterator
      >(stat.n.begin(), stat.id.begin()), 
      tmp_shrt.begin(), // will store the grid cell indices
      tmp_long.begin()  // will store the concentrations per grid cell
    );
    // writing to aux
    mtx::arr<real_t> &n_tot = aux.at("n_tot");
    n_tot(n_tot.ijk) = real_t(0); // as some grid cells could be void of particles
    thrust::counting_iterator<thrust_size_t> iter(0);
    thrust::for_each(iter, iter + (n.first - tmp_shrt.begin()), 
      sdm::copy_from_device<real_t>( 
        grid.nx(), tmp_shrt, tmp_long, n_tot, 
        real_t(1) / grid.dx() / grid.dy() / grid.dz() * si::cubic_metres
      )
    );
  }

  // calculating number of particles with radius greater than a given threshold
  {
    // nested functor
    class threshold_counter
    {
      private: const sdm::stat_t<real_t> &stat;
      private: const thrust_real_t threshold;
      public: threshold_counter(
        const sdm::stat_t<real_t> &stat,
        const thrust_real_t threshold
      ) : stat(stat), threshold(threshold) {}
      public: thrust_size_t operator()(const thrust_size_t id) const
      {
        return stat.xi[id] > threshold ? stat.n[id] : 0;
      }
    };
    // zeroing the temporary var
    thrust::fill(tmp_long.begin(), tmp_long.end(), 0); // TODO: is it needed
    // doing the reduction
    thrust::pair<
      thrust::device_vector<int>::iterator, 
      thrust::device_vector<thrust_size_t>::iterator
    > n = thrust::reduce_by_key(
      stat.ij.begin(), stat.ij.end(),
      thrust::transform_iterator<
        threshold_counter, 
        thrust::device_vector<thrust_size_t>::iterator,
        thrust_size_t
      >(stat.id.begin(), threshold_counter(stat, F_xi->transform(thrust_real_t(500e-9)))),
      tmp_shrt.begin(), // will store the grid cell indices
      tmp_long.begin()  // will store the concentrations per grid cell
    );
    // writing to aux
    mtx::arr<real_t> &n_ccn = aux.at("n_ccn");
    n_ccn(n_ccn.ijk) = real_t(0); // as some grid cells could be void of particles
    thrust::counting_iterator<thrust_size_t> iter(0);
    thrust::for_each(iter, iter + (n.first - tmp_shrt.begin()), 
      sdm::copy_from_device<real_t>( 
        grid.nx(), tmp_shrt, tmp_long, n_ccn, 
        real_t(1) / grid.dx() / grid.dy() / grid.dz() * si::cubic_metres
      )
    );
  }
}

template <typename real_t>
void eqs_todo_sdm<real_t>::sd_advection(
  const quantity<si::time, real_t> dt,
  const mtx::arr<real_t> &Cx,
  const mtx::arr<real_t> &Cy
)
{
  // getting velocities from the Eulerian model (TODO: if constant_velocity -> only once!)
  {
    thrust::counting_iterator<thrust_size_t> iter(0);
    thrust::for_each(iter, iter + envi.vx.size(), 
      sdm::copy_to_device<real_t, thrust_real_t>(
        grid.nx() + 1, Cx, envi.vx, (grid.dx() / si::metres) / (dt / si::seconds)
      )
    );
  }
  {
    thrust::counting_iterator<thrust_size_t> iter(0);
    thrust::for_each(iter, iter + envi.vy.size(), 
      sdm::copy_to_device<real_t, thrust_real_t>(
        grid.nx(), Cy, envi.vy, (grid.dx() / si::metres) / (dt / si::seconds)
      )
    );
  }

  // performing advection using odeint // TODO sedimentation
  F_xy->advance(stat.xy, dt);

  // periodic boundary conditions (in-place transforms)
  thrust::transform(stat.x_begin, stat.x_end, stat.x_begin, 
    sdm::modulo<thrust_real_t>(
      grid.nx() * (grid.dx() / si::metres)
    )
  );
  thrust::transform(stat.y_begin, stat.y_end, stat.y_begin, 
    sdm::modulo<thrust_real_t>(
      grid.ny() * (grid.dy() / si::metres)
    )
  );
}
  
template <typename real_t>
void eqs_todo_sdm<real_t>::sd_condevap(
  const quantity<si::time, real_t> dt
)
{
  // growing/shrinking the droplets
  F_xi->advance(stat.xi, dt);
}
#  endif

// explicit instantiations
#if defined(USE_FLOAT)
template class eqs_todo_sdm<float>;
#endif
#if defined(USE_DOUBLE)
template class eqs_todo_sdm<double>;
#endif
#if defined(USE_LDOUBLE)
template class eqs_todo_sdm<long double>;
#endif
