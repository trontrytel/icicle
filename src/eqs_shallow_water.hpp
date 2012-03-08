/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date January 2012 - February 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 @  @brief contains definition of the eqs_shallow_water class - a system of 2D shallow-water equations
 */
#ifndef EQS_SHALLOW_WATER_HPP
#  define EQS_SHALLOW_WATER_HPP

#  include "cmn.hpp"
#  include "eqs.hpp"
#  include "grd.hpp"

/** @brief the 2D shallow-water equations system
 *
 * Consult chapter 3 in Vallis 2008 for a detailed derivation.
 *
 * The key assumptions are:
 * - horizontal scale is much larger than the vertical scale (\f$ u \approx u(x) \f$)
 * - hydrostatic equillibrium
 * - constant density
 * 
 * Nomenclature:
 * - \f$ \eta(x,y) \f$ - (absolute) height of the fluid surface
 * - \f$ \eta_0(x,y) \f$ - bathymetry
 * - \f$ h = \eta - \eta_0 \f$ - thickness of the fluid layer
 * - \f$ \vec{u} = (u,v) \f$
 * - \f$ \nabla_z = \partial_x + \partial_y \f$
 */
template <typename real_t>
class eqs_shallow_water : public eqs<real_t> 
{
  private: ptr_vector<struct eqs<real_t>::gte> sys;
  public: ptr_vector<struct eqs<real_t>::gte> &system() { return sys; }

  private: struct params
  {
    quantity<si::acceleration, real_t> g;
    quantity<si::length, real_t> h_unit;
    quantity<velocity_times_length, real_t> q_unit;
    int idx_h;
    unique_ptr<mtx::arr<real_t>> dHdx, dHdy;
  };

  // TODO: Coriolis (implicit!)
  // TODO: could the numerics be placed somewhere else...
  // TODO: stencil-extent like method?
  /// @brief Shallow Water Equations: Momentum forcings for the X and Y coordinates
  private: 
  template <int di, int dj>
  class forcings : public rhs<real_t>
  { 
    private: struct params *par;
    private: quantity<si::length, real_t> dxy;
    public: forcings(params &par, quantity<si::length, real_t> dxy) : par(&par), dxy(dxy) {} 
    public: void explicit_part(mtx::arr<real_t> &R, mtx::arr<real_t> **psi, mtx::idx &ijk) 
    { 
      assert((di == 0 && dj == 1) || (di == 1 && dj == 0));

      /// momentum equation:
      /// \f$ \partial_t u + u \cdot \nabla_z u = - \frac{1}{\rho} \nabla_z p \f$ 
      ///
      /// pressure in a column of the constant-density fluid:
      /// \f$ p = p_0 - \rho g z = p_0 + \rho g \cdot (\eta(x) - z) \f$
      ///
      /// mass continuity equation: 
      /// \f$ \partial_t h + \nabla_z (h \cdot u) = 0 \f$
      ///
      /// momentum eq. plus u times mass continuity equation:
      /// \f$ \partial_t (uh)  + \nabla_z (uh) = -g h \nabla_z \eta \f$
      R(ijk) -= 
        par->g * par->h_unit * par->h_unit * si::seconds / par->q_unit / si::metres *
        ((*psi[par->idx_h])(ijk)) *
        (
          (
            ((*psi[par->idx_h])(ijk.i + di, ijk.j + dj, ijk.k)) - 
            ((*psi[par->idx_h])(ijk.i - di, ijk.j - dj, ijk.k))
          ) / (real_t(2) * dxy / si::metres)
          + di * (*par->dHdx)(mtx::idx_ijk(ijk.i, ijk.j, 0)) 
          + dj * (*par->dHdy)(mtx::idx_ijk(ijk.i, ijk.j, 0)) 
        );
    };
  };

  private: params par;
  public: eqs_shallow_water(const grd<real_t> &grid, const ini<real_t> &intcond)
  {
    if (grid.nz() != 1) error_macro("only 1D (X or Y) and 2D (XY) simullations supported")

    par.g = real_t(10.) * si::metres_per_second_squared; // TODO: option (+ constants catalogue!)
    par.h_unit = 1 * si::metres;
    par.q_unit = 1 * si::metres * si::metres / si::seconds;

    // reading topography derivatives from the input file
    mtx::idx_ijk xy(mtx::rng(0, grid.nx()-1), mtx::rng(0, grid.ny()-1), 0); // TODO: not-optimal for MPI (each node has the whole topography!)
    par.dHdx.reset(new mtx::arr<real_t>(xy));
    par.dHdy.reset(new mtx::arr<real_t>(xy));
    intcond.populate_scalar_field("dHdx", par.dHdx->ijk, *(par.dHdx));
    intcond.populate_scalar_field("dHdy", par.dHdy->ijk, *(par.dHdy));

    /// momentum equation for x
    if (grid.nx() != 1)
    {
      sys.push_back(new struct eqs<real_t>::gte({
        "qx", "heigh-integrated specific momentum (x)", 
        this->quan2str(par.q_unit), 
        vector<int>({1, 0, 0})
      }));
      sys.back().rhs_terms.push_back(new forcings<1,0>(par, grid.dx())); 
    }

    /// momentum equation for y
    if (grid.ny() != 1)
    {
      sys.push_back(new struct eqs<real_t>::gte({
        "qy", "heigh-integrated specific momentum (y)", 
        this->quan2str(par.q_unit), 
        vector<int>({0, 1, 0})
      }));
      sys.back().rhs_terms.push_back(new forcings<0,1>(par, grid.dy())); 
    }

    /// The mass continuity equation for a column of height \f$ h = \eta - \eta_0 \f$ of a constant-density fluid:
    /// \f$ \partial_t h + \nabla_z ( \vec{u} h ) = 0 \f$
    sys.push_back(new struct eqs<real_t>::gte({
      "h", "thickness of the fluid layer", 
      this->quan2str(par.h_unit), 
      vector<int>({-1, -1, 0})
    }));
    par.idx_h = sys.size() - 1;
  }
};
#endif
