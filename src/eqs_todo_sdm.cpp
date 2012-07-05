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
#include "eqs_todo_sdm.hpp"

#if defined(USE_BOOST_ODEINT) && defined(USE_THRUST)

#  include "phc_lognormal.hpp" // TODO: not here?
#  include "phc_kappa_koehler.hpp" // TODO: not here?

#  include "sdm_ode_ys.hpp"
#  include <list>
using std::list;

template <typename real_t>
eqs_todo_sdm<real_t>::eqs_todo_sdm(
  const grd<real_t> &grid, 
  const vel<real_t> &velocity,
  map<enum processes, bool> opts,
  enum xi_dfntns xi_dfntn,
  enum ode_algos xy_algo,
  enum ode_algos sd_algo,
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
  envi(grid.nx(), grid.ny()),
  xi_dfntn(xi_dfntn)
{
  // TODO: assert that we use no paralellisation or allow some parallelism!
  // TODO: random seed as an option
  // TODO: option: number of cores to use (via Thrust)

  tmp_shrt.resize(grid.nx() * grid.ny());
  tmp_long.resize(grid.nx() * grid.ny()); // is it used anyehere??? TODO
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

  // auxliary variable for concentration 
  ptr_map_insert(this->aux)("m_0", typename eqs<real_t>::axv({
    "m_0", "<r^0> for r > r_min", this->quan2str(real_t(1)/si::cubic_metres),
    typename eqs<real_t>::invariable(false),
    vector<int>({0, 0, 0})
  }));

  // auxliary variable for mean radius
  ptr_map_insert(this->aux)("m_1", typename eqs<real_t>::axv({
    "m_1", "<r^1> for r > r_min", this->quan2str(si::metres/si::cubic_metres),
    typename eqs<real_t>::invariable(false),
    vector<int>({0, 0, 0})
  }));

  // auxliary variable for mean squared radius
  ptr_map_insert(this->aux)("m_2", typename eqs<real_t>::axv({
    "m_2", "<r^2> for r > r_min", this->quan2str(si::square_metres/si::cubic_metres),
    typename eqs<real_t>::invariable(false),
    vector<int>({0, 0, 0})
  }));

  // auxliary variable for mean cubed radius
  ptr_map_insert(this->aux)("m_3", typename eqs<real_t>::axv({
    "m_3", "<r^3> for r > r_min", this->quan2str(si::cubic_metres/si::cubic_metres),
    typename eqs<real_t>::invariable(false),
    vector<int>({0, 0, 0})
  }));

  // auxliary variable for mean cubed radius
  ptr_map_insert(this->aux)("m_6", typename eqs<real_t>::axv({
    "m_6", "<r^6> for r > r_min", this->quan2str(si::cubic_metres),
    typename eqs<real_t>::invariable(false),
    vector<int>({0, 0, 0})
  }));

  // initialising super-droplets
  sd_init(min_rd, max_rd, mean_rd1, mean_rd2, sdev_rd1, sdev_rd2, n1_tot, n2_tot, sd_conc_mean, kappa);

 typedef odeint::euler<
      thrust::device_vector<real_t>, // state type
      real_t, // value_type
      thrust::device_vector<real_t>, // deriv type
      real_t, // time type
      odeint::thrust_algebra,
      odeint::thrust_operations
    > algo_euler;

  typedef odeint::runge_kutta4<
      thrust::device_vector<real_t>, // state type
      real_t, // value_type
      thrust::device_vector<real_t>, // deriv type
      real_t, // time type
      odeint::thrust_algebra,
      odeint::thrust_operations
    > algo_rk4;

  typedef odeint::modified_midpoint<
      thrust::device_vector<real_t>, // state type
      real_t, // value_type
      thrust::device_vector<real_t>, // deriv type
      real_t, // time type
      odeint::thrust_algebra,
      odeint::thrust_operations
    > algo_mmid;

  // initialising ODE right-hand-sides
  switch (xy_algo) // advection
  {
    case euler: F_xy.reset(new sdm::ode_xy<real_t, algo_euler>(stat, envi)); break;
    case mmid : F_xy.reset(new sdm::ode_xy<real_t, algo_mmid >(stat, envi)); break;
    case rk4  : F_xy.reset(new sdm::ode_xy<real_t, algo_rk4  >(stat, envi)); break;
    default: assert(false);
  }
  switch (xi_algo) // condensation/evaporation
  {
    case euler: switch (xi_dfntn)
    {
      case id : F_xi.reset(new sdm::ode_xi<real_t, algo_euler, sdm::xi_id<real_t>>(stat, envi, grid)); break;
      case ln : F_xi.reset(new sdm::ode_xi<real_t, algo_euler, sdm::xi_ln<real_t>>(stat, envi, grid)); break;
      case p2 : F_xi.reset(new sdm::ode_xi<real_t, algo_euler, sdm::xi_p2<real_t>>(stat, envi, grid)); break;
      case p3 : F_xi.reset(new sdm::ode_xi<real_t, algo_euler, sdm::xi_p3<real_t>>(stat, envi, grid)); break;
      default: assert(false);
    } 
    break;
    case mmid: switch (xi_dfntn)
    {
      case id : F_xi.reset(new sdm::ode_xi<real_t, algo_mmid, sdm::xi_id<real_t>>(stat, envi, grid)); break;
      case ln : F_xi.reset(new sdm::ode_xi<real_t, algo_mmid, sdm::xi_ln<real_t>>(stat, envi, grid)); break;
      case p2 : F_xi.reset(new sdm::ode_xi<real_t, algo_mmid, sdm::xi_p2<real_t>>(stat, envi, grid)); break;
      case p3 : F_xi.reset(new sdm::ode_xi<real_t, algo_mmid, sdm::xi_p3<real_t>>(stat, envi, grid)); break;
      default: assert(false);
    } 
    break;
    case rk4  : switch (xi_dfntn) 
    {
      case id : F_xi.reset(new sdm::ode_xi<real_t, algo_rk4,   sdm::xi_id<real_t>>(stat, envi, grid)); break;
      case ln : F_xi.reset(new sdm::ode_xi<real_t, algo_rk4,   sdm::xi_ln<real_t>>(stat, envi, grid)); break;
      case p2 : F_xi.reset(new sdm::ode_xi<real_t, algo_rk4,   sdm::xi_p2<real_t>>(stat, envi, grid)); break;
      case p3 : F_xi.reset(new sdm::ode_xi<real_t, algo_rk4,   sdm::xi_p3<real_t>>(stat, envi, grid)); break;
      default: assert(false);
    }
    break;
    default: assert(false);
  } 
  switch (sd_algo) // sedimentation
  {
    case euler : switch (xi_dfntn)
    {
      case id : F_ys.reset(new sdm::ode_ys<real_t,algo_euler,sdm::xi_id<real_t>>(stat, envi)); break;
      case ln : F_ys.reset(new sdm::ode_ys<real_t,algo_euler,sdm::xi_ln<real_t>>(stat, envi)); break;
      case p2 : F_ys.reset(new sdm::ode_ys<real_t,algo_euler,sdm::xi_p2<real_t>>(stat, envi)); break;
      case p3 : F_ys.reset(new sdm::ode_ys<real_t,algo_euler,sdm::xi_p3<real_t>>(stat, envi)); break;
      default: assert(false);
    }
    break;
    case mmid: switch (xi_dfntn)
    {
      case id : F_ys.reset(new sdm::ode_ys<real_t,algo_mmid,  sdm::xi_id<real_t>>(stat, envi)); break;
      case ln : F_ys.reset(new sdm::ode_ys<real_t,algo_mmid,  sdm::xi_ln<real_t>>(stat, envi)); break;
      case p2 : F_ys.reset(new sdm::ode_ys<real_t,algo_mmid,  sdm::xi_p2<real_t>>(stat, envi)); break;
      case p3 : F_ys.reset(new sdm::ode_ys<real_t,algo_mmid,  sdm::xi_p3<real_t>>(stat, envi)); break;
      default: assert(false);
    }
    break;
    case rk4   : switch (xi_dfntn)
    {
      case id : F_ys.reset(new sdm::ode_ys<real_t,algo_rk4,  sdm::xi_id<real_t>>(stat, envi)); break;
      case ln : F_ys.reset(new sdm::ode_ys<real_t,algo_rk4,  sdm::xi_ln<real_t>>(stat, envi)); break;
      case p2 : F_ys.reset(new sdm::ode_ys<real_t,algo_rk4,  sdm::xi_p2<real_t>>(stat, envi)); break;
      case p3 : F_ys.reset(new sdm::ode_ys<real_t,algo_rk4,  sdm::xi_p3<real_t>>(stat, envi)); break;
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
  real_t seed = 1234.;

  // initialise particle coordinates
  {
    sdm::rng<real_t> rand_x(0, grid.nx() * (grid.dx() / si::metres), seed);
    thrust::generate(stat.x_begin, stat.x_end, rand_x);
    seed = rand_x();
  }
  {
    sdm::rng<real_t> rand_y(0, grid.ny() * (grid.dy() / si::metres), seed);
    thrust::generate(stat.y_begin, stat.y_end, rand_y);
    seed = rand_y();
  }

  // initialise particle dry size spectrum 
  // TODO: assert that the distribution is < epsilon at rd_min and rd_max
  thrust::generate( // rd3 temporarily means logarithms of radius!
    stat.rd3.begin(), stat.rd3.end(), 
    sdm::rng<real_t>(log(min_rd), log(max_rd), seed) 
  ); 
  {
    real_t multi = log(max_rd/min_rd) / sd_conc_mean 
      * (grid.dx() * grid.dy() * grid.dz() / si::cubic_metres); 
    thrust::transform(
      stat.rd3.begin(), stat.rd3.end(), 
      stat.n.begin(), 
      sdm::lognormal<real_t>(
        real_t(mean_rd1) * si::metres, 
        real_t(sdev_rd1), 
        real_t(multi * n1_tot) / si::cubic_metres, 
        real_t(mean_rd2) * si::metres, 
        real_t(sdev_rd2), 
        real_t(multi * n2_tot) / si::cubic_metres
      ) 
    );
    // converting rd back from logarihms to rd3 
    struct exp3x { real_t operator()(real_t x) { return exp(3*x); } };
    thrust::transform(stat.rd3.begin(), stat.rd3.end(), stat.rd3.begin(), exp3x());
  }

  // initialise kappas
  thrust::fill(stat.kpa.begin(), stat.kpa.end(), kappa);
}

// sorting out which particle belongs to which grid cell
template <typename real_t>
void eqs_todo_sdm<real_t>::sd_ij()
{
  // computing stat.i 
  thrust::transform(
    stat.x_begin, stat.x_end,
    stat.i.begin(),
    sdm::divide_by_constant<real_t>(grid.dx() / si::metres)
  );
  // computing stat.j
  thrust::transform(
    stat.y_begin, stat.y_end,
    stat.j.begin(),
    sdm::divide_by_constant<real_t>(grid.dy() / si::metres)
  );
  // computing stat.ij
  thrust::transform(
    stat.i.begin(), stat.i.end(),
    stat.j.begin(),
    stat.ij.begin(),
    sdm::ravel_indices(grid.ny())
  );
}

template <typename real_t>
void eqs_todo_sdm<real_t>::sd_sync_in(
  const mtx::arr<real_t> &rhod,
  const mtx::arr<real_t> &rhod_th,
  const mtx::arr<real_t> &rhod_rv
)
{
  // getting thermodynamic fields from the Eulerian model 
  thrust::counting_iterator<thrust_size_t> iter(0);
  {
    thrust::for_each(iter, iter + envi.rhod.size(), 
      sdm::blitz2thrust<real_t, real_t>(
        grid.nx(), rhod, envi.rhod
      )
    ); // TODO: this could be done just once in the kinematic model!
  }
  {
    thrust::for_each(iter, iter + envi.rhod_th.size(), 
      sdm::blitz2thrust<real_t, real_t>(
        grid.nx(), rhod_th, envi.rhod_th 
      )
    );
  }
  {
    thrust::for_each(iter, iter + envi.rhod_rv.size(), 
      sdm::blitz2thrust<real_t, real_t>(
        grid.nx(), rhod_rv, envi.rhod_rv
      )
    );
  }

  // calculating the derived fields (T,p,r)
  // TODO! should be done just once!!! TODO - now it's repeated in ode_xi::adjust()
  class rpT
  {
    private: sdm::envi_t<real_t> &envi;
    public: rpT(sdm::envi_t<real_t> &envi) : envi(envi) {}
    public: void operator()(const thrust_size_t ij) 
    { 
      envi.r[ij] = envi.rhod_rv[ij] / envi.rhod[ij]; 
      envi.p[ij] = phc::p<real_t>(
        envi.rhod_th[ij] * si::kelvins * si::kilograms / si::cubic_metres, 
        quantity<phc::mixing_ratio, real_t>(envi.r[ij])
      ) / si::pascals; 
      envi.T[ij] = phc::T<real_t>(
        (envi.rhod_th[ij] / envi.rhod[ij]) * si::kelvins, 
        envi.p[ij] * si::pascals, 
        quantity<phc::mixing_ratio, real_t>(envi.r[ij])
      ) / si::kelvins; 
    }
  };
  thrust::for_each(iter, iter + envi.n_cell, rpT(envi));
}

template <typename real_t>
void eqs_todo_sdm<real_t>::sd_sync_out(
  mtx::arr<real_t> &rhod_th,
  mtx::arr<real_t> &rhod_rv
)
{
  thrust::counting_iterator<thrust_size_t> zero(0);
  thrust::sequence(tmp_shrt.begin(), tmp_shrt.end()); // TODO: could be optimised...
  thrust::for_each(zero, zero + envi.n_cell, 
    sdm::thrust2blitz<real_t>(grid.nx(), tmp_shrt, envi.rhod_th, rhod_th)
  );
  thrust::for_each(zero, zero + envi.n_cell,
    sdm::thrust2blitz<real_t>(grid.nx(), tmp_shrt, envi.rhod_rv, rhod_rv)
  );
}

// computing diagnostics
template <typename real_t>
void eqs_todo_sdm<real_t>::sd_diag(ptr_unordered_map<string, mtx::arr<real_t>> &aux)
{
  // making a copy of ij
  thrust::copy(stat.ij.begin(), stat.ij.end(), stat.sorted_ij.begin());
  // filling-in stat.sorted_id
  thrust::sequence(stat.sorted_id.begin(), stat.sorted_id.end());
  // sorting stat.ij and stat.sorted_id
  thrust::sort_by_key(
    stat.sorted_ij.begin(), stat.sorted_ij.end(),
    stat.sorted_id.begin()
  );

  // calculating super-droplet concentration (per grid cell)
  {
    // zeroing the temporary var 
    thrust::fill(tmp_long.begin(), tmp_long.end(), 0); // TODO: is it needed?
    // doing the reduction
    thrust::pair<
      thrust::device_vector<int>::iterator, 
      thrust::device_vector<thrust_size_t>::iterator
    > n = thrust::reduce_by_key(
      stat.sorted_ij.begin(), stat.sorted_ij.end(),
      thrust::make_constant_iterator(1),
      tmp_shrt.begin(), // will store the grid cell indices
      tmp_long.begin()  // will store the concentrations per grid cell
    );
    // writing to aux
    mtx::arr<real_t> &sd_conc = aux.at("sd_conc");
    sd_conc(sd_conc.ijk) = real_t(0); // as some grid cells could be void of particles
    thrust::counting_iterator<thrust_size_t> iter(0);
    thrust::for_each(iter, iter + (n.first - tmp_shrt.begin()), 
      sdm::thrust2blitz<real_t>(
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
      stat.sorted_ij.begin(), stat.sorted_ij.end(),
      thrust::permutation_iterator<
        thrust::device_vector<thrust_size_t>::iterator,
        thrust::device_vector<thrust_size_t>::iterator
      >(stat.n.begin(), stat.sorted_id.begin()), 
      tmp_shrt.begin(), // will store the grid cell indices
      tmp_long.begin()  // will store the concentrations per grid cell
    );
    // writing to aux
    mtx::arr<real_t> &n_tot = aux.at("n_tot");
    n_tot(n_tot.ijk) = real_t(0); // as some grid cells could be void of particles
    thrust::counting_iterator<thrust_size_t> iter(0);
    thrust::for_each(iter, iter + (n.first - tmp_shrt.begin()), 
      sdm::thrust2blitz<real_t>( 
        grid.nx(), tmp_shrt, tmp_long, n_tot, 
        real_t(1) / grid.dx() / grid.dy() / grid.dz() * si::cubic_metres
      )
    );
  }

  // calculating k-th moment of the particle distribution (for particles with radius greater than a given threshold)
  for (int &k : list<int>({0,1,2,3,6}))
  {
    // zeroing the temporary var
    thrust::fill(tmp_real.begin(), tmp_real.end(), 0); // TODO: is it needed

    // doing the reduction, i.e. 
    thrust::pair<
      decltype(tmp_shrt.begin()),
      decltype(tmp_real.begin()) 
    > n;

    switch (xi_dfntn)
    {
      case p2:
      n = thrust::reduce_by_key(
        stat.sorted_ij.begin(), stat.sorted_ij.end(),
        thrust::transform_iterator<
          sdm::moment_counter<real_t, sdm::xi_p2<real_t>>,
          decltype(stat.sorted_id.begin()),
          real_t
        >(stat.sorted_id.begin(), sdm::moment_counter<real_t, sdm::xi_p2<real_t>>(stat, real_t(500e-9), k)),//
        tmp_shrt.begin(), // will store the grid cell indices
        tmp_real.begin()  // will store the concentrations per grid cell
      );
      break;
      default: assert(false); //TODO: ln, p3, id, ...
    }

    // writing to aux
    ostringstream tmp;
    tmp << "m_" << k;
    mtx::arr<real_t> &out = aux.at(tmp.str());
    out(out.ijk) = real_t(0); // as some grid cells could be void of particles
    thrust::counting_iterator<thrust_size_t> iter(0);
    thrust::for_each(iter, iter + (n.first - tmp_shrt.begin()), 
      sdm::thrust2blitz<real_t>( 
        grid.nx(), tmp_shrt, tmp_real, out, 
        real_t(1) / grid.dx() / grid.dy() / grid.dz() * si::cubic_metres
      )
    );
  }

  // 
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
      sdm::blitz2thrust<real_t, real_t>(
        grid.nx() + 1, Cx, envi.vx, (grid.dx() / si::metres) / (dt / si::seconds)
      )
    );
  }
  {
    thrust::counting_iterator<thrust_size_t> iter(0);
    thrust::for_each(iter, iter + envi.vy.size(), 
      sdm::blitz2thrust<real_t, real_t>(
        grid.nx(), Cy, envi.vy, (grid.dx() / si::metres) / (dt / si::seconds)
      )
    );
  }

  // performing advection using odeint 
  F_xy->advance(stat.xy, dt);
}

template <typename real_t>
void eqs_todo_sdm<real_t>::sd_sedimentation(
  const quantity<si::time, real_t> dt,
  const mtx::arr<real_t> &rhod
)
{
  // performing advection using odeint
  F_ys->advance(stat.xy, dt);
}

template <typename real_t>
void eqs_todo_sdm<real_t>::sd_periodic_x()
{
  // periodic boundary condition (in-place transforms)
  thrust::transform(stat.x_begin, stat.x_end, stat.x_begin, 
    sdm::modulo<real_t>(
      grid.nx() * (grid.dx() / si::metres)
    )
  );
}

template <typename real_t>
void eqs_todo_sdm<real_t>::sd_periodic_y()
{
  // periodic boundary condition (in-place transforms)
  thrust::transform(stat.y_begin, stat.y_end, stat.y_begin, 
    sdm::modulo<real_t>(
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
  F_xi->advance(stat.xi, dt, 20); // TODO: maximal timestep as an option!
  assert(*thrust::min_element(stat.xi.begin(), stat.xi.end()) > 0);
}
#endif

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
