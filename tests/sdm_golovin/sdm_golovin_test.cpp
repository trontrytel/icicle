/** @file
 *  @example sdm_golovin/sdm_golovin_test.cpp
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date October 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @section DESCRIPTION
 *    Tests the implementation of the Monte-Carlo collision scheme 
 *    see fig. 2 in Shima et al. 2009 
 */

#include "../../src/sdm/sdm.hpp"
#include "../../src/grd_carthesian.hpp"
#include "../../src/vel_func_uniform.hpp"
#include <list>
using std::list;

typedef double real_t;

int main()
{
  // setting from Shima et al. 2009 (sec. 5.1.4)
  quantity<si::length, real_t> 
    dx = 1e3 * si::metres,
    dy = 1e3 * si::metres,
    dz = 1 * si::metres; 

  // TODO: try with different nx, ny
  int nx, ny, nz;
  nx = ny = nz = 1;

  //
  grd_carthesian<real_t> grid(dx, dy, dz, nx, ny, nz);
  vel_func_uniform<real_t> velocity(grid, 
    0 * si::metres_per_second,
    0 * si::metres_per_second,
    0 * si::metres_per_second
  );
 
  map<enum sdm::processes, bool> opts({
    {sdm::adve, false},
    {sdm::cond, false},
    {sdm::sedi, false},
    {sdm::coal, true },
    {sdm::chem, false},
  });

  int sd_conc_mean = 256; 
  real_t min_rd = .1e-6, max_rd = 100e-6;

  // TODO: to be changed into an exponential distro
  real_t 
    mean_rd1 = .04e-6,
    mean_rd2 = .15e-6,
    sdev_rd1 = 1.4,
    sdev_rd2 = 1.6,
    n1_tot = 60e6,
    n2_tot = 40e6,
    kappa = 0.61;

  map<enum sdm::chem_gas, quantity<phc::mixing_ratio, real_t>> opt_gas({
    {sdm::gSO2,  0},
    {sdm::gO3,   0},
    {sdm::gH2O2, 0},
  });
  map<enum sdm::chem_aq, quantity<si::mass, real_t>> opt_aq({
    {sdm::H,    0 * si::kilograms},
    {sdm::OH,   0 * si::kilograms},
    {sdm::SO2,  0 * si::kilograms},
    {sdm::O3,   0 * si::kilograms},
    {sdm::H2O2, 0 * si::kilograms},
    {sdm::HSO3, 0 * si::kilograms},
    {sdm::SO3,  0 * si::kilograms}
  });

/*

  // TODO: for(openmp, cuda)
/*
  sdm::sdm<real_t, sdm::openmp> particles(
    grid,
    velocity,
    opts,
    sdm::p2,
    sdm::euler, sdm::euler, sdm::euler, sdm::euler,
    1, 1, 1, 1, 1, // substeps
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
  );
*/
}
