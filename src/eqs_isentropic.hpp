/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date February 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 @  @brief contains definition of the eqs_isentropic class - a system of 3D isentropic equations
 */
#ifndef EQS_ISENTROPIC_HPP
#  define EQS_ISENTROPIC_HPP

#  include "cmn.hpp"
#  include "eqs.hpp"
#  include "grd.hpp"

/// @brief the 3D isentropic equations system
template <typename real_t>
class eqs_isentropic : public eqs<real_t> 
{
  private: ptr_vector<struct eqs<real_t>::gte> sys;
  public: ptr_vector<struct eqs<real_t>::gte> &system() { return sys; }

  private: struct params
  {
    quantity<si::acceleration, real_t> g;
    quantity<specific_heat_capacity, real_t> cp;
    quantity<si::dimensionless, real_t> Rd_over_cp;
    quantity<si::pressure, real_t> p_unit, p0;
    quantity<velocity_times_pressure, real_t> q_unit;
    int nlev;
    vector<int> idx_dp;
    vector<quantity<si::temperature, real_t>> theta;
  };

  private: 
  template <int di, int dj>
  class forcings : public rhs<real_t>
  { 
    private: struct params *par;
    private: quantity<si::length, real_t> dxy;
    private: int lev;
    public: forcings(params &par, quantity<si::length, real_t> dxy, int lev) 
      : par(&par), dxy(dxy), lev(lev) {} 
    public: void operator()(mtx::arr<real_t> &R, mtx::arr<real_t> **psi, mtx::idx &ijk) 
    { 
      assert((di == 0 && dj > 0) || (di > 0 && dj == 0));

      // we need a local "halo" for the M spatial derivative
      mtx::rng 
        ii(ijk.lbound(mtx::i) - 1, ijk.ubound(mtx::i) + 1), 
        jj(ijk.lbound(mtx::j) - 1, ijk.ubound(mtx::j) + 1),
        ll(0, par->nlev);

      // allocating temporary arrays TODO: should be initialised only once per simulation
      mtx::arr<real_t> M(mtx::idx_ijk(ii, jj, ijk.k)); // mid-layer Montgomery potentials 
      mtx::arr<real_t> p(mtx::idx_ijk(ii, jj, ll)); // pressures at layer boundries

      // calculating pressure (TODO: should be done once per timestep)
      p(mtx::idx_ijk(ii, jj, ll)) = real_t(0); // TODO: parameter
      for (int l = par->nlev - 1; l >= 0; --l)
        p(mtx::idx_ijk(ii, jj, l)) = 
          (*psi[par->idx_dp[l]])(mtx::idx_ijk(ii, jj, ijk.k)) 
          + p(mtx::idx_ijk(ii, jj, l+1));

      // "specific" Montgomery potential
      quantity<specific_energy, real_t> M_unit = 1 * si::metres_per_second_squared * si::metres;

      // calculating Montgomery potential (TODO: should be done incrementaly)
      M(mtx::idx_ijk(ii,jj,ijk.k)) = real_t(0); // TODO: topography
      for (int l = 0; l < par->nlev; ++l)
        M(mtx::idx_ijk(ii,jj,ijk.k)) += par->theta[lev] * par->cp * pow(par->p_unit / par->p0, par->Rd_over_cp) 
           / M_unit // to make it dimensionless
           * pow( real_t(.5) * (p(mtx::idx_ijk(ii, jj, l)) + p(mtx::idx_ijk(ii, jj, l+1))), par->Rd_over_cp );

      // eq. (6) in Szmelter & Smolarkiewicz 2011, Computers & Fluids
      R(ijk) -= 
        M_unit * par->p_unit / (real_t(2 * (di + dj)) * dxy) // real units of the rhs
        * si::seconds / par->q_unit // to get a dimensionless forcing
        * ((*psi[par->idx_dp[lev]])(ijk)) * (
          (M(ijk.i + di, ijk.j + dj, ijk.k))
          - 
          (M(ijk.i - di, ijk.j - dj, ijk.k))
        );
    };
  };

  private: params par;
  public: eqs_isentropic(const grd<real_t> &grid, int nlev, quantity<si::acceleration, real_t> g)
  {
    if (grid.nz() != 1) error_macro("TODO") // TODO!

    // constants
    par.g = g;
    par.p_unit = 1 * si::pascals;
    par.q_unit = 1 * si::pascals * si::metres / si::seconds;
    par.nlev = nlev;
    par.p0 = 100000 * si::pascals; // definition of potential temperature
    par.cp = 1005 * si::joules / si::kilograms / si::kelvins;
    par.Rd_over_cp = si::constants::codata::R / (0.02896 * si::kilograms / si::moles) / par.cp; // TODO!!!!

    // potential temperature profile
    par.theta.resize(nlev);
    par.theta[0] = 200 * si::kelvins; // TEMP
    par.theta[1] = 300 * si::kelvins; // TEMP
    par.theta[2] = 400 * si::kelvins; // TEMP

    // dp var indices (to be filled in in the loop below)
    par.idx_dp.resize(nlev);

    for (int lint = 0; lint < nlev; ++lint)
    {
      string lstr = boost::lexical_cast<std::string>(lint);

      sys.push_back(new struct eqs<real_t>::gte({
        "dp_" + lstr, "pressure thickness of fluid layer " + lstr,
        this->quan2str(par.p_unit), 
        vector<int>({-1, -1, 0}),
        typename eqs<real_t>::groupid(lint)
      }));
      par.idx_dp[lint] = sys.size() - 1;

      sys.push_back(new struct eqs<real_t>::gte({
        "qx_" + lstr, "layer-integrated specific momentum (x) of fluid layer " + lstr, 
        this->quan2str(par.q_unit), 
        vector<int>({1, 0, 0}),
        typename eqs<real_t>::groupid(lint)
      }));
      sys.back().source_terms.push_back(new forcings<1,0>(par, grid.dx(), lint)); 

      sys.push_back(new struct eqs<real_t>::gte({
        "qy_" + lstr, "layer-integrated specific momentum (y) of fluid layer " + lstr, 
        this->quan2str(par.q_unit), 
        vector<int>({0, 1, 0}),
        typename eqs<real_t>::groupid(lint)
      }));
      sys.back().source_terms.push_back(new forcings<0,1>(par, grid.dy(), lint)); 
    }
  }
};
#endif
