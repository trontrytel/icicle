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
    quantity<si::pressure, real_t> p_unit, p0, p_max;
    quantity<velocity_times_pressure, real_t> q_unit;
    vector<int> idx_dp;
    unique_ptr<mtx::arr<real_t>> dtheta, dHdx, dHdy, M, p; // TODO: M and p need halo_filling - should be handled somehow by the solver!
  };

  private: 
  template <int di, int dj>
  class forcings : public rhs<real_t>
  { 
    private: struct params *par;
    private: quantity<si::length, real_t> dxy;
    private: int lev;
    private: bool calc_p, calc_M;
    public: forcings(params &par, quantity<si::length, real_t> dxy, int lev, bool calc_p, bool calc_M) 
      : par(&par), dxy(dxy), lev(lev), calc_p(calc_p), calc_M(calc_M) {} 
    public: void operator()(mtx::arr<real_t> &R, mtx::arr<real_t> **psi, mtx::idx &ijk) 
    { 
      assert((di == 0 && dj == 1) || (di == 1 && dj == 0));
      int nlev = par->idx_dp.size();

      // we need a local "halo" for the M spatial derivative
      mtx::rng 
        ii(ijk.lbound(mtx::i) - 1, ijk.ubound(mtx::i) + 1), 
        jj(ijk.lbound(mtx::j) - 1, ijk.ubound(mtx::j) + 1),
        ll(0, nlev);

      // calculating pressures at the surfaces (done once per timestep) // TODO: need fill_halos!!!
      if (calc_p) // TODO: will not work in parallel!!!
      {
cerr << "calculating pressures!" << endl;
        (*par->p)(mtx::idx_ijk(ii, jj, ll)) = par->p_max / si::pascals;
        for (int l = nlev - 1; l >= 0; --l)
          (*par->p)(mtx::idx_ijk(ii, jj, l)) = 
            (*psi[par->idx_dp[l]])(mtx::idx_ijk(ii, jj, ijk.k)) 
            + (*par->p)(mtx::idx_ijk(ii, jj, l+1));
      }

      // "specific" Montgomery potential
      quantity<specific_energy, real_t> M_unit = 1 * si::metres_per_second_squared * si::metres;

      // calculating Montgomery potential (TODO: should be done incrementaly)
      if (calc_p) 
        (*par->M)(mtx::idx_ijk(ii,jj,ijk.k)) = real_t(0); // TODO: parallel!!!
      if (calc_M) 
      {
cerr << "calculating Montgomery at lev=" << lev << endl;
        // this assumes we go from the bottom up to the top layer
        // (the order of equations in the ctor below)
        (*par->M)(mtx::idx_ijk(ii,jj,ijk.k)) += par->cp * pow(par->p_unit / par->p0, par->Rd_over_cp) 
           / M_unit * si::kelvins * // to make it dimensionless
           ((*par->dtheta)(mtx::idx_ijk(lev)))(0,0,0) * // TODO: nicer syntax needed!
           real_t(.5) * // Exner function within a layer as an average of surface Ex-fun values
           ( 
             pow((*par->p)(mtx::idx_ijk(ii, jj, lev  )), par->Rd_over_cp) +
             pow((*par->p)(mtx::idx_ijk(ii, jj, lev+1)), par->Rd_over_cp)
           );
      }

      // eq. (6) in Szmelter & Smolarkiewicz 2011, Computers & Fluids
      R(ijk) -= 
        M_unit * par->p_unit / si::metres * // real units of the rhs
        si::seconds / par->q_unit * // inv. unit of R (to get a dimensionless forcing)
        ((*psi[par->idx_dp[lev]])(ijk)) * 
        ( 
          par->g / si::metres_per_second_squared *
          ( // di = 0 || dj = 0
            di * (*par->dHdx)(mtx::idx_ijk(ijk.i, ijk.j, 0)) +
            dj * (*par->dHdy)(mtx::idx_ijk(ijk.i, ijk.j, 0))
          ) 
          +
          ( // spatial derivative
            ((*par->M)(ijk.i + di, ijk.j + dj, ijk.k)) - 
            ((*par->M)(ijk.i - di, ijk.j - dj, ijk.k))
          ) / (real_t(2) * dxy / si::metres)
        );
    };
  };

  private: params par;
  public: eqs_isentropic(const grd<real_t> &grid, const ini<real_t> &intcond,
    int nlev, 
    quantity<si::pressure, real_t> p_max,
    quantity<si::acceleration, real_t> g
  )
  {
    if (grid.nz() != 1) error_macro("TODO") // TODO!

    // constants
    par.g = g;
    par.p_max = p_max;
    par.p_unit = 1 * si::pascals;
    par.q_unit = 1 * si::pascals * si::metres / si::seconds;
    par.p0 = 100000 * si::pascals; // definition of potential temperature
    par.cp = 1005 * si::joules / si::kilograms / si::kelvins;
    par.Rd_over_cp = si::constants::codata::R / (0.02896 * si::kilograms / si::moles) / par.cp; // TODO!!!!

    // potential temperature profile
    par.dtheta.reset(new mtx::arr<real_t>(mtx::idx_ijk(mtx::rng(0,nlev-1))));
    intcond.populate_scalar_field("dtheta", par.dtheta->ijk, *(par.dtheta));

    // topography
    mtx::idx_ijk xy(mtx::rng(0, grid.nx()-1), mtx::rng(0, grid.ny()-1), 0); // TODO: not-optimal for MPI (each node has the whole topography!)
    par.dHdx.reset(new mtx::arr<real_t>(xy));
    par.dHdy.reset(new mtx::arr<real_t>(xy));
    intcond.populate_scalar_field("dHdx", par.dHdx->ijk, *(par.dHdx));
    intcond.populate_scalar_field("dHdy", par.dHdy->ijk, *(par.dHdy));

    // allocating temporary arrays (only once per simulation)
    mtx::rng // TODO: non-optimal for MPI!
      ii(0 - 1, grid.nx() + 1),
      jj(0 - 1, grid.ny() + 1),
      ll(0, nlev); // this means nlev+1 surfaces
    par.M.reset(new mtx::arr<real_t>(mtx::idx_ijk(ii, jj, 0))); // mid-layer Montgomery potentials
    par.p.reset(new mtx::arr<real_t>(mtx::idx_ijk(ii, jj, ll))); // pressures at layer boundries       

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
      sys.back().source_terms.push_back(new forcings<1,0>(par, grid.dx(), lint, lint==0, true)); 

      sys.push_back(new struct eqs<real_t>::gte({
        "qy_" + lstr, "layer-integrated specific momentum (y) of fluid layer " + lstr, 
        this->quan2str(par.q_unit), 
        vector<int>({0, 1, 0}),
        typename eqs<real_t>::groupid(lint)
      }));
      sys.back().source_terms.push_back(new forcings<0,1>(par, grid.dy(), lint, false, false)); 
    }
  }
};
#endif
