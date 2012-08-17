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
#  include "phc_terminal_vel.hpp" // TODO: not here?
#  include "phc_kappa_koehler.hpp" // TODO: not here?

#  include "sdm_ode_cond.hpp"
#  include "sdm_ode_adve.hpp"
#  include "sdm_ode_sedi.hpp"
#  include "sdm_ode_chem.hpp"
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
  map<enum sdm::chem_gas, quantity<phc::mixing_ratio, real_t>> opt_gas,
  map<enum sdm::chem_aq, quantity<si::mass, real_t>> opt_aq
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

  // temporary space with nx * ny size
  tmp_shrt.resize(grid.nx() * grid.ny());
  tmp_long.resize(grid.nx() * grid.ny()); 
  tmp_real.resize(grid.nx() * grid.ny());

  // temporary space with n_part size
  tmp_long_id.resize(stat.n_part);
  tmp_real_id.resize(stat.n_part);

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

  // auxliary variable for mean SO3 density within a droplet
  ptr_map_insert(this->aux)("c_SO3", typename eqs<real_t>::axv({
    "c_SO3", "<c_SO3> for r > r_min", this->quan2str(si::kilograms / si::cubic_metres),
    typename eqs<real_t>::invariable(false),
    vector<int>({0, 0, 0})
  }));

  // initialising super-droplets
  sd_init(min_rd, max_rd, mean_rd1, mean_rd2, sdev_rd1, sdev_rd2, n1_tot, n2_tot, sd_conc_mean, kappa, 
          opt_gas, opt_aq);

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

  // TODO: adams_bashforth_moulton...

  // initialising ODE right-hand-sides
  switch (xy_algo) // advection
  {
    case euler: F_adve.reset(new sdm::ode_adve<real_t, algo_euler>(stat, envi)); break;
    case mmid : F_adve.reset(new sdm::ode_adve<real_t, algo_mmid >(stat, envi)); break;
    case rk4  : F_adve.reset(new sdm::ode_adve<real_t, algo_rk4  >(stat, envi)); break;
    default: assert(false);
  }
  switch (xi_algo) // condensation/evaporation
  {
    case euler: switch (xi_dfntn)
    {
      case id : F_cond.reset(new sdm::ode_cond<real_t, algo_euler, sdm::xi_id<real_t>>(stat, envi, grid)); break;
      case ln : F_cond.reset(new sdm::ode_cond<real_t, algo_euler, sdm::xi_ln<real_t>>(stat, envi, grid)); break;
      case p2 : F_cond.reset(new sdm::ode_cond<real_t, algo_euler, sdm::xi_p2<real_t>>(stat, envi, grid)); break;
      case p3 : F_cond.reset(new sdm::ode_cond<real_t, algo_euler, sdm::xi_p3<real_t>>(stat, envi, grid)); break;
      default: assert(false);
    } 
    break;
    case mmid: switch (xi_dfntn)
    {
      case id : F_cond.reset(new sdm::ode_cond<real_t, algo_mmid, sdm::xi_id<real_t>>(stat, envi, grid)); break;
      case ln : F_cond.reset(new sdm::ode_cond<real_t, algo_mmid, sdm::xi_ln<real_t>>(stat, envi, grid)); break;
      case p2 : F_cond.reset(new sdm::ode_cond<real_t, algo_mmid, sdm::xi_p2<real_t>>(stat, envi, grid)); break;
      case p3 : F_cond.reset(new sdm::ode_cond<real_t, algo_mmid, sdm::xi_p3<real_t>>(stat, envi, grid)); break;
      default: assert(false);
    } 
    break;
    case rk4  : switch (xi_dfntn) 
    {
      case id : F_cond.reset(new sdm::ode_cond<real_t, algo_rk4,   sdm::xi_id<real_t>>(stat, envi, grid)); break;
      case ln : F_cond.reset(new sdm::ode_cond<real_t, algo_rk4,   sdm::xi_ln<real_t>>(stat, envi, grid)); break;
      case p2 : F_cond.reset(new sdm::ode_cond<real_t, algo_rk4,   sdm::xi_p2<real_t>>(stat, envi, grid)); break;
      case p3 : F_cond.reset(new sdm::ode_cond<real_t, algo_rk4,   sdm::xi_p3<real_t>>(stat, envi, grid)); break;
      default: assert(false);
    }
    break;
    default: assert(false);
  } 
  switch (sd_algo) // sedimentation
  {
    case euler : switch (xi_dfntn)
    {
      case id : F_sedi.reset(new sdm::ode_sedi<real_t,algo_euler,sdm::xi_id<real_t>>(stat, envi)); break;
      case ln : F_sedi.reset(new sdm::ode_sedi<real_t,algo_euler,sdm::xi_ln<real_t>>(stat, envi)); break;
      case p2 : F_sedi.reset(new sdm::ode_sedi<real_t,algo_euler,sdm::xi_p2<real_t>>(stat, envi)); break;
      case p3 : F_sedi.reset(new sdm::ode_sedi<real_t,algo_euler,sdm::xi_p3<real_t>>(stat, envi)); break;
      default: assert(false);
    }
    break;
    case mmid: switch (xi_dfntn)
    {
      case id : F_sedi.reset(new sdm::ode_sedi<real_t,algo_mmid,  sdm::xi_id<real_t>>(stat, envi)); break;
      case ln : F_sedi.reset(new sdm::ode_sedi<real_t,algo_mmid,  sdm::xi_ln<real_t>>(stat, envi)); break;
      case p2 : F_sedi.reset(new sdm::ode_sedi<real_t,algo_mmid,  sdm::xi_p2<real_t>>(stat, envi)); break;
      case p3 : F_sedi.reset(new sdm::ode_sedi<real_t,algo_mmid,  sdm::xi_p3<real_t>>(stat, envi)); break;
      default: assert(false);
    }
    break;
    case rk4   : switch (xi_dfntn)
    {
      case id : F_sedi.reset(new sdm::ode_sedi<real_t,algo_rk4,  sdm::xi_id<real_t>>(stat, envi)); break;
      case ln : F_sedi.reset(new sdm::ode_sedi<real_t,algo_rk4,  sdm::xi_ln<real_t>>(stat, envi)); break;
      case p2 : F_sedi.reset(new sdm::ode_sedi<real_t,algo_rk4,  sdm::xi_p2<real_t>>(stat, envi)); break;
      case p3 : F_sedi.reset(new sdm::ode_sedi<real_t,algo_rk4,  sdm::xi_p3<real_t>>(stat, envi)); break;
      default: assert(false);
    }
    break;
    default: assert(false);
  }
  switch (chem_algo) // chemistry
  {
    case euler : switch (xi_dfntn)
    {
      case id : F_chem.reset(new sdm::ode_chem<real_t,algo_euler,sdm::xi_id<real_t>>(stat, envi, opt_gas, opt_aq)); break;
      case ln : F_chem.reset(new sdm::ode_chem<real_t,algo_euler,sdm::xi_ln<real_t>>(stat, envi, opt_gas, opt_aq)); break;
      case p2 : F_chem.reset(new sdm::ode_chem<real_t,algo_euler,sdm::xi_p2<real_t>>(stat, envi, opt_gas, opt_aq)); break;
      case p3 : F_chem.reset(new sdm::ode_chem<real_t,algo_euler,sdm::xi_p3<real_t>>(stat, envi, opt_gas, opt_aq)); break;
      default: assert(false);
    }
    break;
    case mmid: switch (xi_dfntn)
    {
      case id : F_chem.reset(new sdm::ode_chem<real_t,algo_mmid,  sdm::xi_id<real_t>>(stat, envi, opt_gas, opt_aq)); break;
      case ln : F_chem.reset(new sdm::ode_chem<real_t,algo_mmid,  sdm::xi_ln<real_t>>(stat, envi, opt_gas, opt_aq)); break;
      case p2 : F_chem.reset(new sdm::ode_chem<real_t,algo_mmid,  sdm::xi_p2<real_t>>(stat, envi, opt_gas, opt_aq)); break;
      case p3 : F_chem.reset(new sdm::ode_chem<real_t,algo_mmid,  sdm::xi_p3<real_t>>(stat, envi, opt_gas, opt_aq)); break;
      default: assert(false);
    }
    break;
    case rk4   : switch (xi_dfntn)
    {
      case id : F_chem.reset(new sdm::ode_chem<real_t,algo_rk4,  sdm::xi_id<real_t>>(stat, envi, opt_gas, opt_aq)); break;
      case ln : F_chem.reset(new sdm::ode_chem<real_t,algo_rk4,  sdm::xi_ln<real_t>>(stat, envi, opt_gas, opt_aq)); break;
      case p2 : F_chem.reset(new sdm::ode_chem<real_t,algo_rk4,  sdm::xi_p2<real_t>>(stat, envi, opt_gas, opt_aq)); break;
      case p3 : F_chem.reset(new sdm::ode_chem<real_t,algo_rk4,  sdm::xi_p3<real_t>>(stat, envi, opt_gas, opt_aq)); break;
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
  const real_t kappa,
  map<enum sdm::chem_gas, quantity<phc::mixing_ratio, real_t>> opt_gas,
  map<enum sdm::chem_aq, quantity<si::mass, real_t>> opt_aq
)
{
  // initialise particle coordinates
  {
    sdm::uniform_random_real<real_t> rand_x(0, grid.nx() * (grid.dx() / si::metres), seed);
    thrust::generate(stat.x_begin, stat.x_end, rand_x);
    seed = rand_x();
  }
  {
    sdm::uniform_random_real<real_t> rand_y(0, grid.ny() * (grid.dy() / si::metres), seed);
    thrust::generate(stat.y_begin, stat.y_end, rand_y);
    seed = rand_y();
  }

  // initialise particle dry size spectrum 
  // TODO: assert that the distribution is < epsilon at rd_min and rd_max
  thrust::generate( // rd3 temporarily means logarithms of radius!
    stat.rd3.begin(), stat.rd3.end(), 
    sdm::uniform_random_real<real_t>(log(min_rd), log(max_rd), seed) 
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

  // initialise chemical components
  thrust::fill(
    stat.c_aq.begin() +  sdm::H      * stat.n_part, 
    stat.c_aq.begin() + (sdm::H + 1) * stat.n_part, 
    opt_aq[sdm::H] / si::kilograms
  );
  thrust::fill(
    stat.c_aq.begin() +  sdm::OH      * stat.n_part, 
    stat.c_aq.begin() + (sdm::OH + 1) * stat.n_part, 
    opt_aq[sdm::OH] / si::kilograms
  );
  thrust::fill(
    stat.c_aq.begin() +  sdm::SO2      * stat.n_part, 
    stat.c_aq.begin() + (sdm::SO2 + 1) * stat.n_part, 
    opt_aq[sdm::SO2] / si::kilograms
  );
  thrust::fill(
    stat.c_aq.begin() +  sdm::O3      * stat.n_part,  
    stat.c_aq.begin() + (sdm::O3 + 1) * stat.n_part, 
    opt_aq[sdm::O3] / si::kilograms
  );
  thrust::fill(
    stat.c_aq.begin() +  sdm::H2O2      * stat.n_part, 
    stat.c_aq.begin() + (sdm::H2O2 + 1) * stat.n_part, 
    opt_aq[sdm::H2O2] / si::kilograms
  );
  thrust::fill(
    stat.c_aq.begin() +  sdm::HSO3      * stat.n_part,
    stat.c_aq.begin() + (sdm::HSO3 + 1) * stat.n_part,
    opt_aq[sdm::HSO3] / si::kilograms
  );
  thrust::fill(
    stat.c_aq.begin() +  sdm::SO3      * stat.n_part,
    stat.c_aq.begin() + (sdm::SO3 + 1) * stat.n_part, 
    opt_aq[sdm::SO3]  / si::kilograms
  );
  thrust::fill(
    stat.c_aq.begin() +  sdm::HSO4      * stat.n_part,
    stat.c_aq.begin() + (sdm::HSO4 + 1) * stat.n_part, 
    0. //TODO?
  );
  thrust::fill(
    stat.c_aq.begin() +  sdm::SO4      * stat.n_part,
    stat.c_aq.begin() + (sdm::SO4 + 1) * stat.n_part, 
    0.//TODO?
  );
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

template <typename real_t>
void eqs_todo_sdm<real_t>::sd_sort()
{
  // filling-in stat.sorted_id with a sequence
  thrust::sequence(stat.sorted_id.begin(), stat.sorted_id.end());
  // making a copy of ij
  thrust::copy(stat.ij.begin(), stat.ij.end(), stat.sorted_ij.begin());
  // sorting stat.ij and stat.sorted_id
  thrust::sort_by_key(
    stat.sorted_ij.begin(), stat.sorted_ij.end(),
    stat.sorted_id.begin()
  );
}

template <typename real_t>
void eqs_todo_sdm<real_t>::sd_shuffle_and_sort()
{
  // filling-in sorted_id with a random sequence 
  // that's a try to get the std::random_shuffle functionality using parallelisable Thrust calls
  {
    // filling-in stat.sorted_id with a sequence (TODO: is it neccesarry here?)
    thrust::sequence(stat.sorted_id.begin(), stat.sorted_id.end());
    // using sorted_ij as temporary space for the random key
    thrust::device_vector<int> &random_key = stat.sorted_ij; 
    {
      sdm::uniform_random_int<int> rand(INT_MIN, INT_MAX, seed); // TODO: c++ version of INT_MIN/INT_MAX
      thrust::generate(random_key.begin(), random_key.end(), rand);
      seed = rand();
    }
    thrust::stable_sort_by_key( // sorting with random key!
      random_key.begin(), random_key.end(),
      stat.sorted_id.begin()
    );
  }
  // making a permutated copy of ij
  thrust::permutation_iterator<
    decltype(stat.ij.begin()), // values
    decltype(stat.sorted_id.begin()) // indices
  > perm(stat.ij.begin(), stat.sorted_id.begin());
  thrust::copy(perm, perm + stat.n_part, stat.sorted_ij.begin()); 
  // sorting stat.ij and stat.sorted_id
  thrust::stable_sort_by_key( // TODO is stable needed?
    stat.sorted_ij.begin(), stat.sorted_ij.end(),
    stat.sorted_id.begin()
  );
}

// computing diagnostics
template <typename real_t>
void eqs_todo_sdm<real_t>::sd_diag(ptr_unordered_map<string, mtx::arr<real_t>> &aux)
{
  // calculating super-droplet concentration (per grid cell)
  {
    // TODO: repeated in sd_coalescence!
    // zeroing the temporary var 
    thrust::fill(tmp_long.begin(), tmp_long.end(), 0); // TODO: is it needed?
    // doing the reduction
    thrust::pair<
      decltype(tmp_shrt.begin()), 
      decltype(tmp_long.begin())
    > n = thrust::reduce_by_key(
      stat.sorted_ij.begin(), stat.sorted_ij.end(),
      thrust::make_constant_iterator(1),
      tmp_shrt.begin(), // will store the grid cell indices
      tmp_long.begin()  // will store the concentrations per grid cell
    );
    // writing to aux
    mtx::arr<real_t> &sd_conc = aux.at("sd_conc");
    sd_conc(sd_conc.ijk) = real_t(0); // as some grid cells could be void of particles
    thrust::counting_iterator<thrust_size_t> zero(0);
    thrust::for_each(zero, zero + (n.first - tmp_shrt.begin()), 
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
        >(stat.sorted_id.begin(), sdm::moment_counter<real_t, sdm::xi_p2<real_t>>(stat, real_t(500e-9), k)),// TODO:!!! 500e-9 as an option!!!
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

  // calculating chemical stuff
  if (true/* TODO opt[chem]*/)
  {
    typedef pair<enum sdm::chem_aq, string> keyval;
    for (keyval &kv : list<keyval>({{sdm::SO3,"c_SO3"}}))
    {
      // zeroing the temporary var
      thrust::fill(tmp_real.begin(), tmp_real.end(), 0); // TODO: is it needed

      // doing the reduction, i.e. 
      thrust::pair<
        decltype(tmp_shrt.begin()),
        decltype(tmp_real.begin())
      > n;

      n = thrust::reduce_by_key(
        stat.sorted_ij.begin(), stat.sorted_ij.end(),
        thrust::transform_iterator<
          sdm::chem_counter<real_t>,
          decltype(stat.sorted_id.begin()),
          real_t
        >(stat.sorted_id.begin(), sdm::chem_counter<real_t>(stat, kv.first)),
        tmp_shrt.begin(), // will store the grid cell indices
        tmp_real.begin()  // will store the concentrations per grid cell
      );

      // writing to aux
      mtx::arr<real_t> &out = aux.at(kv.second);
      out(out.ijk) = real_t(0); // as some grid cells could be void of particles
      thrust::counting_iterator<thrust_size_t> iter(0);
      thrust::for_each(iter, iter + (n.first - tmp_shrt.begin()),
        sdm::thrust2blitz<real_t>(
          grid.nx(), tmp_shrt, tmp_real, out,
          real_t(1) / grid.dx() / grid.dy() / grid.dz() * si::cubic_metres
        )
      );
    }
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
  F_adve->advance(stat.xy, dt);
}

template <typename real_t>
void eqs_todo_sdm<real_t>::sd_sedimentation(
  const quantity<si::time, real_t> dt,
  const mtx::arr<real_t> &rhod
)
{
  // performing advection using odeint
  F_sedi->advance(stat.xy, dt);
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
void eqs_todo_sdm<real_t>::sd_recycle() 
{
  // invalidate out-of-domain droplets
  struct invalidate
  {
    // member fields
    eqs_todo_sdm<real_t> &nest;
    const real_t Y;
    // ctor
    invalidate(eqs_todo_sdm<real_t> &nest) :
      nest(nest), Y(nest.grid.nx() * (nest.grid.dx() / si::metres))
    {}
    // overloaded operator invoked by Thrust
    void operator()(thrust_size_t id)
    {
      if (
        nest.stat.xy[nest.stat.n_part + id] < 0 || 
        nest.stat.xy[nest.stat.n_part + id] > Y
      ) 
      {
        nest.stat.xy[nest.stat.n_part + id] = Y/2.; // TODO
        nest.stat.n[id] = 0; // invalidating
      }
    }
  };
  thrust::counting_iterator<thrust_size_t> zero(0);
  thrust::for_each(zero, zero + stat.n_part, invalidate(*this));

  // TODO: recycle droplets!
}
  
template <typename real_t>
void eqs_todo_sdm<real_t>::sd_condevap(
  const quantity<si::time, real_t> dt
)
{
  // growing/shrinking the droplets
  F_cond->advance(stat.xi, dt, 50); // TODO: maximal timestep as an option!
  assert(*thrust::min_element(stat.xi.begin(), stat.xi.end()) > 0);
}

template <typename real_t>
void eqs_todo_sdm<real_t>::sd_chem(
  const quantity<si::time, real_t> dt
)
{
  F_chem->advance(stat.c_aq, dt, 10); // TODO !!!!!
  assert(*thrust::min_element(stat.c_aq.begin(), stat.c_aq.end()) >= 0);
}

template <typename real_t>
struct eqs_todo_sdm<real_t>::detail
{
  // 
  template <class xi>
  struct collider : xi
  { 
    // member fields
    eqs_todo_sdm<real_t> &nest;
    const quantity<si::time, real_t> dt;
    const quantity<si::volume, real_t> dv;
    sdm::uniform_random_real<real_t> rand;

    // ctor
    collider(
      eqs_todo_sdm<real_t> &nest,
      const quantity<si::time, real_t> dt
    ) : 
      nest(nest), 
      dt(dt),
      dv(nest.grid.dx() * nest.grid.dy() * nest.grid.dz()),
      rand(0, 1, nest.seed)
    {}
   
    // dtor
    ~collider()
    {
      nest.seed = rand();
    }

    // overloaded operator invoked by Thrust
    void operator()(sdm::thrust_size_t key)
    {
      // we take into account every second droplet
      if (key % 2 != 0) return; 

      // and we only work with droplets within the same cell
      if (nest.stat.sorted_ij[key] != nest.stat.sorted_ij[key+1]) return;

      // we get the id-s and ij-s from the sorted random permutation
      sdm::thrust_size_t
        id1 = nest.stat.sorted_id[key],
        id2 = nest.stat.sorted_id[key+1],
        ij = nest.stat.sorted_ij[key];

      // chosing id1 to correspond to the higher multiplicity
      if (nest.stat.n[id1] < nest.stat.n[id2]) std::swap(id1, id2);

      // evaluating the probability
      real_t prob = 
        real_t(nest.stat.n[id1]) // max(n1,n2) multiplicity (targe/projectile configuration)
        * real_t(nest.tmp_real[ij]) // scaling factor
        * abs(
          phc::vt<real_t>(
            this->rw_of_xi(nest.stat.xi[id1]) * si::metres,
            nest.envi.T[ij] * si::kelvins,
            nest.envi.rhod[ij] * si::kilograms / si::cubic_metres // that's the dry air density
          ) -
          phc::vt<real_t>(
            this->rw_of_xi(nest.stat.xi[id2]) * si::metres,
            nest.envi.T[ij] * si::kelvins,
            nest.envi.rhod[ij] * si::kilograms / si::cubic_metres // that's the dry air density
          )
        ) // velocity difference
        * dt // timestep
        / dv // volume
        * si::square_metres * phc::pi<real_t>() * pow(this->rw_of_xi(nest.stat.xi[id1]) + this->rw_of_xi(nest.stat.xi[id2]), 2) // geometrical cross-section
        * real_t(10); // collection efficiency TODO!!!
      assert(prob < 1);

      // tossing a random number and returning if unlucky
      if (prob < rand()) return;

      // performing the coalescence event if lucky
      if (nest.stat.n[id1] != nest.stat.n[id2])
      {
        // multiplicity change (eq. 12 in Shima et al. 2009)
        nest.stat.n[id1] -= nest.stat.n[id2]; 
        // wet radius change (eq. 13 in Shima et al. 2009)
        nest.stat.xi[id2] = this->xi_of_rw3(
          this->rw3_of_xi(nest.stat.xi[id1]) + 
          this->rw3_of_xi(nest.stat.xi[id2])
        ); 
        // dry radius change (eq. 13 in Shima et al. 2009)
        nest.stat.rd3[id2] = nest.stat.rd3[id1] + nest.stat.rd3[id2];
      }
      else
      {
        // TODO: eqs. 16-19 in Shima et al. 2009
        assert(false);
      }
    }
  };
};

template <typename real_t>
void eqs_todo_sdm<real_t>::sd_coalescence(
  const quantity<si::time, real_t> dt
)
{
  thrust::counting_iterator<thrust_size_t> zero(0);

  // first getting the number of super-droplets per cell, 
  // a list of cells to consider, and probability scale-factors
  {
    // zeroing the temporary var 
    thrust::fill(tmp_long.begin(), tmp_long.end(), 0); // TODO: is it needed?
    // doing the reduction
    thrust::pair<
      decltype(tmp_shrt.begin()), 
      decltype(tmp_long.begin())
    > n = thrust::reduce_by_key(
      stat.sorted_ij.begin(), stat.sorted_ij.end(),
      thrust::make_constant_iterator(1),
      tmp_shrt.begin(), // will store the grid cell indices
      tmp_long.begin()  // will store the concentrations per grid cell
    ); 
  
    // filling tmp_real with (n*(n-1)/2)/floor(n/2) for given ij
    thrust::fill(tmp_real.begin(), tmp_real.end(), real_t(0));

    struct scale_factor  
    {
      eqs_todo_sdm<real_t> &nest;
      scale_factor(eqs_todo_sdm<real_t> &nest) : nest(nest) {}
      inline real_t sf(thrust_size_t n) { return real_t((n*(n-1))/2) / (n/2); }
      void operator()(thrust_size_t i)
      {
        nest.tmp_real[nest.tmp_shrt[i]] = sf(nest.tmp_long[i]);
      }
    };
    thrust::for_each(zero, zero + (n.first - tmp_shrt.begin()), scale_factor(*this));
  }

  // getting particle number in a sequence through a segmented prefix sum 
  thrust::exclusive_scan_by_key(
    stat.sorted_ij.begin(), stat.sorted_ij.end(),
    stat.sorted_id.begin(),
    tmp_long_id.begin() 
  );

  switch (xi_dfntn)
  {
    case p2:
      // the collider uses n and n+1, hence stopping at n_part - 1
      thrust::for_each(zero, zero + stat.n_part - 1, typename detail::collider<sdm::xi_p2<real_t>>(*this, dt));
      break;
    default: assert(false); //TODO: ln, p3, id, ...
  }
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
