/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date April-June 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @brief definition of the eqs_todo_sdm class
 */

#include "eqs_todo_sdm.hpp"
#if defined(USE_THRUST) && defined(USE_BOOST_ODEINT)
#  include "../sdm/sdm.hpp"

template <typename real_t>
struct eqs_todo_sdm<real_t>::detail
{
  // both sdm<openmp> and sdm<cuda> inherit from sdm_proto 
  unique_ptr<sdm::sdm_proto<real_t>> particles;
};

template <typename real_t>
eqs_todo_sdm<real_t>::eqs_todo_sdm(
  const grd<real_t> &grid, 
  const vel<real_t> &velocity,
  map<enum sdm::processes, bool> opts,
  const enum sdm::xi_dfntns xi_dfntn,
  const enum sdm::ode_algos adve_algo,
  const enum sdm::ode_algos sedi_algo,
  const enum sdm::ode_algos cond_algo,
  const enum sdm::ode_algos chem_algo,
  const int adve_sstp,
  const int sedi_sstp,
  const int cond_sstp,
  const int chem_sstp,
  const int coal_sstp,
  const real_t sd_conc_mean,
  const real_t min_rd,
  const real_t max_rd, 
  const real_t mean_rd1, // dry aerosol initial distribution parameters
  const real_t mean_rd2,
  const real_t sdev_rd1,
  const real_t sdev_rd2,
  const real_t n1_tot,
  const real_t n2_tot,
  const real_t kappa,
  map<enum sdm::chem_gas, quantity<phc::mixing_ratio, real_t>> opt_gas,
  map<enum sdm::chem_aq, quantity<si::mass, real_t>> opt_aq
) : eqs_todo<real_t>(grid, &this->par), pimpl(new detail())
{
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

  // auxliary variable for mean SO2 density within a droplet
  ptr_map_insert(this->aux)("c_SO3", typename eqs<real_t>::axv({
    "c_SO3", "<c_SO3> for r > r_min", this->quan2str(si::kilograms / si::cubic_metres),
    typename eqs<real_t>::invariable(false),
    vector<int>({0, 0, 0})
  }));

  // initialising super-droplets
  if (true) // TODO!
    pimpl->particles.reset(new sdm::sdm<real_t, sdm::openmp>(
      grid,
      velocity,
      opts,
      xi_dfntn,
      adve_algo, sedi_algo, cond_algo, chem_algo,
      adve_sstp, sedi_sstp, cond_sstp, chem_sstp, coal_sstp,
      sd_conc_mean,
      min_rd,
      max_rd, 
      mean_rd1, // dry aerosol initial distribution parameters
      mean_rd2, // TODO: encapsulate in a map/enum
      sdev_rd1,
      sdev_rd2,
      n1_tot,
      n2_tot,
      kappa,
      opt_gas, 
      opt_aq 
));
/*
  else
    pimpl->particles.reset(new sdm::sdm<real_t, sdm::cuda>(
      grid, 
      velocity,
      opts,
      xi_dfntn,
      adve_algo, sedi_algo, cond_algo, chem_algo,
      adve_sstp, sedi_sstp, cond_sstp, chem_sstp, coal_sstp,
      sd_conc_mean,
      min_rd,
      max_rd, 
      mean_rd1, // dry aerosol initial distribution parameters
      mean_rd2, // TODO: encapsulate in a map/enum
      sdev_rd1,
      sdev_rd2,
      n1_tot,
      n2_tot,
      kappa,
      opt_gas, 
      opt_aq 
));
*/
}

template <typename real_t>
void eqs_todo_sdm<real_t>::adjustments(
  int n, // TODO: mo¿e jednak bez n...
  vector<ptr_vector<mtx::arr<real_t>>> &psi,
  ptr_unordered_map<string, mtx::arr<real_t>> &aux, 
  const ptr_vector<mtx::arr<real_t>> C,
  const quantity<si::time, real_t> dt,
  bool record
)
{
  pimpl->particles->adjustments(
    psi[this->par.idx_rhod_th][n],
    psi[this->par.idx_rhod_rv][n],
    aux, C, dt, record
  );
}
#endif

// explicit instantiations
#define ICICLE_INSTANTIATE_CLASS eqs_todo_sdm
#include "../cmn/cmn_instant.hpp"
