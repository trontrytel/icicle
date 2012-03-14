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
  private: struct params
  {
    quantity<si::acceleration, real_t> g;
    quantity<si::length, real_t> h_unit;
    quantity<velocity_times_length, real_t> q_unit;
    quantity<si::dimensionless, real_t> dHdxy_unit;
    int idx_h;
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
    private: int idx_dHdxy;

    // ctor
    public: forcings(params &par, quantity<si::length, real_t> dxy, int idx_dHdxy) : 
      par(&par), dxy(dxy), idx_dHdxy(idx_dHdxy) 
    {
      assert((di == 0 && dj == 1) || (di == 1 && dj == 0));
    } 

    public: void explicit_part(
      mtx::arr<real_t> &R, 
      mtx::arr<real_t> **aux,
      mtx::arr<real_t> **psi
    ) 
    { 
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
      R(R.ijk) -= 
        par->g * par->h_unit * par->h_unit * si::seconds / par->q_unit / si::metres *
        ((*psi[par->idx_h])(R.ijk)) *
        (
          (
            ((*psi[par->idx_h])(mtx::idx_ijk(R.ijk.i + di, R.ijk.j + dj, 0))) - 
            ((*psi[par->idx_h])(mtx::idx_ijk(R.ijk.i - di, R.ijk.j - dj, 0)))
          ) / (real_t(2) * dxy / si::metres)
          + 
          (*aux[idx_dHdxy])(mtx::idx_ijk(R.ijk.i, R.ijk.j, 0)) 
        );
    };
  };

  private: params par;

  private: ptr_vector<struct eqs<real_t>::gte> sys;
  public: ptr_vector<struct eqs<real_t>::gte> &system() { return sys; }

  private: ptr_vector<struct eqs<real_t>::axv> aux;
  public: ptr_vector<struct eqs<real_t>::axv> &auxvars() { return aux; }

  // ctor
  public: eqs_shallow_water(const grd<real_t> &grid)
  {
    if (grid.nz() != 1) error_macro("only 1D (X or Y) and 2D (XY) simullations supported")

    par.g = real_t(10.) * si::metres_per_second_squared; // TODO: option (+ constants catalogue!)
    par.h_unit = 1 * si::metres;
    par.q_unit = 1 * si::metres * si::metres / si::seconds;
    par.dHdxy_unit = 1;

    // topography derivatives
    int idx_dHdx = -1, idx_dHdy = -1; 
    {   
      if (grid.nx() != 1)
      {   
        aux.push_back(new struct eqs<real_t>::axv({
          "dHdx", "spatial derivative of the topography (X)", this->quan2str(par.dHdxy_unit),
          vector<int>({0, 0, 1}) // dimspan
        }));
        idx_dHdx = aux.size() - 1;
      }   
      if (grid.ny() != 1)
      {   
        aux.push_back(new struct eqs<real_t>::axv({
          "dHdy", "spatial derivative of the topograpy (Y)", this->quan2str(par.dHdxy_unit),
          vector<int>({0, 0, 1}) // dimspan
        }));
        idx_dHdy = aux.size() - 1;
      }
    }

    /// momentum equation for x
    if (grid.nx() != 1)
    {
      sys.push_back(new struct eqs<real_t>::gte({
        "qx", "heigh-integrated specific momentum (x)", 
        this->quan2str(par.q_unit), 
        typename eqs<real_t>::positive_definite(false),
        vector<int>({1, 0, 0})
      }));
      sys.back().rhs_terms.push_back(new forcings<1,0>(par, grid.dx(), idx_dHdx)); 
    }

    /// momentum equation for y
    if (grid.ny() != 1)
    {
      sys.push_back(new struct eqs<real_t>::gte({
        "qy", "heigh-integrated specific momentum (y)", 
        this->quan2str(par.q_unit), 
        typename eqs<real_t>::positive_definite(false),
        vector<int>({0, 1, 0})
      }));
      sys.back().rhs_terms.push_back(new forcings<0,1>(par, grid.dy(), idx_dHdy)); 
    }

    /// The mass continuity equation for a column of height \f$ h = \eta - \eta_0 \f$ of a constant-density fluid:
    /// \f$ \partial_t h + \nabla_z ( \vec{u} h ) = 0 \f$
    sys.push_back(new struct eqs<real_t>::gte({
      "h", "thickness of the fluid layer", 
      this->quan2str(par.h_unit), 
      typename eqs<real_t>::positive_definite(true),
      vector<int>({-1, -1, 0})
    }));
    par.idx_h = sys.size() - 1;
  }
};
#endif
