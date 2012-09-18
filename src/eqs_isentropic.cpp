/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date February 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 @  @brief contains definition of the eqs_isentropic class - a system of 3D isentropic equations
 */

#include "rhs_explicit.hpp"
#include "eqs_isentropic.hpp"
#include "cmn/cmn_error.hpp"

#include <boost/lexical_cast.hpp>

#include <boost/assign/ptr_map_inserter.hpp>
using boost::assign::ptr_map_insert;

template <typename real_t>
class eqs_isentropic<real_t>::rayleigh_damping : public rhs_implicit<real_t> // (aka sponge layer, aka gravity-wave absorber)
{
  private: real_t tau(int lev, int nlev, int abslev, real_t absamp) const
  {
    if (lev <= abslev) return real_t(0);
    real_t pi = acos(real_t(-1));
    return absamp * pow(sin(.5*pi * (lev - abslev) / ((nlev - 1) - abslev)), 2);
  }
  private: real_t C;
  public: rayleigh_damping(int lev, int nlev, int abslev, real_t absamp)
    : C(-tau(lev, nlev, abslev, absamp))
  { }
  public: real_t implicit_part(const quantity<si::time, real_t>) const
  { return C; }
};

template <typename real_t>
template <int di, int dj>
class eqs_isentropic<real_t>::montgomery_grad : public rhs_explicit<real_t>
{ 
  // members
  private: struct params *par;
  private: quantity<si::length, real_t> dxy;
  private: int lev, nlev;
  private: string idx_dHdxy;
  private: bool calc_p, calc_M;

  // ctor
  public: montgomery_grad(params &par, quantity<si::length, real_t> dxy, 
    const string &idx_dHdxy, int lev, bool calc_p, bool calc_M
  ) :
    par(&par), dxy(dxy), idx_dHdxy(idx_dHdxy), lev(lev), 
    nlev(par.idx_dp.size()), calc_p(calc_p), calc_M(calc_M) 
  {
    assert((di == 0 && dj == 1) || (di == 1 && dj == 0));
  } 
 
  // method
  public: void explicit_part(
    mtx::arr<real_t> &R, 
    ptr_unordered_map<string, mtx::arr<real_t>> &aux,
    const mtx::arr<real_t> * const * const psi,
    const quantity<si::time, real_t> 
  ) const
  {
    const mtx::arr<real_t> &p = aux.at("p"); 
    const mtx::arr<real_t> &M = aux.at("M");
    const mtx::rng &ii = p.ijk.i, &jj = p.ijk.j;

    // calculating pressures at the surfaces (done once per timestep)
    if (calc_p)
    {
      p(mtx::idx_ijk(ii, jj, nlev)) = par->p_top / si::pascals;
      for (int l = nlev - 1; l >= 0; --l)
        p(mtx::idx_ijk(ii, jj, l)) = 
          (*psi[par->idx_dp[l]])(mtx::idx_ijk(ii, jj, 0)) 
          + 
          p(mtx::idx_ijk(ii, jj, l+1));
      assert(isfinite(sum(p)));
    }

    // calculating Montgomery potential
    // (assuming we go from the bottom up to the top layer -
    //  the order of equations is defined in the eqs_isentropic ctor below)
    if (calc_M) 
    {
      if (calc_p) 
        M(M.ijk) = phc::c_pd<real_t>() * real_t(phc::exner<real_t>(par->p_unit))
          / par->M_unit // to make it dimensionless
          * par->theta_frst 
          * real_t(.5) 
          * ( 
             pow(p(mtx::idx_ijk(ii, jj, 0)), phc::R_d_over_c_pd<real_t>()) +
             pow(p(mtx::idx_ijk(ii, jj, 1)), phc::R_d_over_c_pd<real_t>())
          );
      else 
      {
        const mtx::arr<real_t> &dtheta(aux.at("dtheta"));
        M(M.ijk) += phc::c_pd<real_t>() * real_t(phc::exner<real_t>(par->p_unit))
          / par->M_unit * si::kelvins * // to make it dimensionless
          (dtheta(mtx::idx_ijk(lev - 1)))(0,0,0) * // TODO: nicer syntax needed!
           real_t(.5) * // Exner function within a layer as an average of surface Ex-fun values
           ( 
             pow(p(mtx::idx_ijk(ii, jj, lev  )), phc::R_d_over_c_pd<real_t>()) +
             pow(p(mtx::idx_ijk(ii, jj, lev+1)), phc::R_d_over_c_pd<real_t>())
           );
      }
      assert(isfinite(sum(M)));
    }

    // eq. (6) in Szmelter & Smolarkiewicz 2011, Computers & Fluids
    const mtx::arr<real_t> &dHdxy(aux.at(idx_dHdxy));
    R(R.ijk) -= 
      par->M_unit * par->p_unit / si::metres * // real units of the rhs
      si::seconds / par->q_unit * // inv. unit of R (to get a dimensionless forcing)
      ((*psi[par->idx_dp[lev]])(R.ijk)) * 
      ( 
        phc::g<real_t>() / si::metres_per_second_squared *
        dHdxy(mtx::idx_ijk(R.ijk.i, R.ijk.j, 0)) 
        +
        ( // spatial derivative
          (M(mtx::idx_ijk(R.ijk.i + di, R.ijk.j + dj, 0))) - 
          (M(mtx::idx_ijk(R.ijk.i - di, R.ijk.j - dj, 0)))
        ) / (real_t(2) * dxy / si::metres)
      );
  };
};

template <typename real_t>
eqs_isentropic<real_t>::eqs_isentropic(const grd<real_t> &grid,
  int nlev, 
  quantity<si::pressure, real_t> p_top,
  quantity<si::temperature, real_t> theta_frst,
  int abslev, real_t absamp
)
{
  if (grid.nz() != 1) error_macro("only 1D (X or Y) or 2D (XY) simulations supported") 

  // parameters
  par.p_top = p_top;
  par.theta_frst = theta_frst;
 
  // units
  par.p_unit = 1 * si::pascals;
  par.q_unit = 1 * si::pascals * si::metres / si::seconds;
  par.theta_unit = 1 * si::kelvins;
  par.M_unit = 1 * si::metres_per_second_squared * si::metres; // "specific" Montgomery potential
  par.dHdxy_unit = 1;

  // potential temperature profile
  ptr_map_insert(this->aux)("dtheta", typename eqs<real_t>::axv({
    "dtheta", "mid-layer potential temperature increment", this->quan2str(par.theta_unit),
    typename eqs<real_t>::invariable(true),
    vector<int>({nlev - 1, 1, 1})
  }));

  // topography
  {
    if (grid.nx() != 1)
    {
      ptr_map_insert(this->aux)("dHdx", typename eqs<real_t>::axv({
        "dHdx", "spatial derivative of the topography (X)", this->quan2str(par.dHdxy_unit),
        typename eqs<real_t>::invariable(true),
        vector<int>({0, 0, 1}) // dimspan
      }));
    }
    if (grid.ny() != 1)
    {
      ptr_map_insert(this->aux)("dHdy", typename eqs<real_t>::axv({
        "dHdy", "spatial derivative of the topograpy (Y)", this->quan2str(par.dHdxy_unit),
        typename eqs<real_t>::invariable(true),
        vector<int>({0, 0, 1}) // dimspan
      }));
    }
  }

  // pressure 
  ptr_map_insert(this->aux)("p", typename eqs<real_t>::axv({
    "p", "pressure", this->quan2str(par.p_unit),
    typename eqs<real_t>::invariable(false),
    vector<int>({0, 0, nlev+1 }), // dimspan
    vector<int>({1, 1, 0}) // halo extent
  }));

  // the Montgomery potential
  ptr_map_insert(this->aux)("M", typename eqs<real_t>::axv({
    "M", "Montgomery potential", this->quan2str(par.M_unit),
    typename eqs<real_t>::invariable(false),
    vector<int>({0, 0, 1}), // dimspan
    vector<int>({1, 1, 0}) // halo extent
  }));
    
  // dp var indices (to be filled in in the loop below)
  par.idx_dp.resize(nlev);

  // loop over fluid layers
  for (int lint = 0; lint < nlev; ++lint)
  {
    string lstr = boost::lexical_cast<std::string>(lint);

    this->sys.push_back(new struct eqs<real_t>::gte({
      "dp_" + lstr, "pressure thickness of fluid layer " + lstr,
      this->quan2str(par.p_unit), 
      eqs<real_t>::positive_definite(true),
      vector<int>({-1, -1, 0}),
      typename eqs<real_t>::groupid(lint)
    }));
    par.idx_dp[lint] = this->sys.size() - 1;

    if (grid.nx() != 1)
    {
      this->sys.push_back(new struct eqs<real_t>::gte({
        "qx_" + lstr, "layer-integrated specific momentum (x) of fluid layer " + lstr, 
        this->quan2str(par.q_unit), 
        eqs<real_t>::positive_definite(false),
        vector<int>({1, 0, 0}),
        typename eqs<real_t>::groupid(lint)
      }));
      this->sys.back().rhs_terms.push_back(
        new montgomery_grad<1,0>(par, grid.dx(), "dHdx", lint, lint==0, true)
      ); 
      if (lint > abslev)
      {
        this->sys.back().rhs_terms.push_back(
          new rayleigh_damping(lint, nlev, abslev, absamp)
        ); 
      }
    }

    if (grid.ny() != 1)
    {
      this->sys.push_back(new struct eqs<real_t>::gte({
        "qy_" + lstr, "layer-integrated specific momentum (y) of fluid layer " + lstr, 
        this->quan2str(par.q_unit), 
        eqs<real_t>::positive_definite(false),
        vector<int>({0, 1, 0}),
        typename eqs<real_t>::groupid(lint)
      }));
      this->sys.back().rhs_terms.push_back(
        new montgomery_grad<0,1>(par, grid.dy(), "dHdy", lint, false, false)
      ); 
      if (lint > abslev)
      {
        this->sys.back().rhs_terms.push_back(
          new rayleigh_damping(lint, nlev, abslev, absamp)
        ); 
      }
    }
  }
  this->init_maps();
}

// explicit instantiations
#define ICICLE_INSTANTIATE_CLASS eqs_isentropic
#include "cmn/cmn_instant.hpp"
