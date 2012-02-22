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

/// @brief the 2D shallow-water equations system
template <typename real_t>
class eqs_shallow_water : public eqs<real_t> 
{
  private: ptr_vector<struct eqs<real_t>::gte> sys;
  public: ptr_vector<struct eqs<real_t>::gte> &system() { return sys; }

  private: struct params
  {
    quantity<si::acceleration, real_t> g;
    quantity<si::length, real_t> dx, dy, h_unit;
    quantity<velocity_times_length, real_t> q_unit;
    int idx_h;
  };

  // TODO: Coriolis
  // TODO: the numerics should be placed somewhere else...
  /// @brief Shallow Water Equations: Momentum forcings for the X coordinate
  private: class forcings_qx : public rhs<real_t>
  { 
    private: struct params *par;
    public: forcings_qx(params &par) : par(&par) {} 
    public: void operator()(mtx::arr<real_t> &R, mtx::arr<real_t> **psi, mtx::idx &ijk) 
    { 
      R(ijk) -= 
        real_t(par->g / (real_t(2) * par->dx) * par->h_unit * par->h_unit * si::seconds / par->q_unit) 
        * ((*psi[par->idx_h])(ijk)) * 
        (
          ((*psi[par->idx_h])(ijk.i + 1, ijk.j, ijk.k) /* TODO: + H0 */)
          - 
          ((*psi[par->idx_h])(ijk.i - 1, ijk.j, ijk.k) /* TODO: + H0 */)
        );
    };
  };

  /// @brief Shallow Water Equations: Momentum forcings for the Y coordinate
  private: class forcings_qy : public rhs<real_t>
  { 
    private: struct params *par;
    public: forcings_qy(params &par) : par(&par) {} 
    public: void operator()(mtx::arr<real_t> &R, mtx::arr<real_t> **psi, mtx::idx &ijk) 
    { 
      R(ijk) -= real_t(par->g / (real_t(2) * par->dy) * par->h_unit * par->h_unit * si::seconds / par->q_unit) 
        * ((*psi[par->idx_h])(ijk)) * 
        (
          ((*psi[par->idx_h])(ijk.i, ijk.j + 1, ijk.k) /* TODO: + H0 */)
          - 
          ((*psi[par->idx_h])(ijk.i, ijk.j - 1, ijk.k) /* TODO: + H0 */)
        );
    };
  };

  private: quantity<si::length, real_t> h_unit;
  private: quantity<velocity_times_length, real_t> q_unit;
  private: params par;
  public: eqs_shallow_water(const grd<real_t> &grid)
  {
    if (grid.nz() != 1) error_macro("only 1D (X or Y) and 2D (XY) simullations supported")

    par.g = real_t(10.) * si::metres_per_second_squared, grid.dx(); // TODO: option
    par.dx = grid.dx();
    par.dy = grid.dy();
    par.h_unit = 1 * si::metres;
    par.q_unit = 1 * si::metres * si::metres / si::seconds;

    sys.push_back(new struct eqs<real_t>::gte({
      "qx", "heigh-integrated specific momentum (x)", 
      this->quan2str(q_unit), 
      vector<int>({1, 0, 0})
    }));
    sys.back().source_terms.push_back(new forcings_qx(par)); 

    sys.push_back(new struct eqs<real_t>::gte({
      "qy", "heigh-integrated specific momentum (y)", 
      this->quan2str(q_unit), 
      vector<int>({0, 1, 0})
    }));
    sys.back().source_terms.push_back(new forcings_qy(par)); 

    sys.push_back(new struct eqs<real_t>::gte({
      "h", "thickness of the fluid layer", 
      this->quan2str(h_unit), 
      vector<int>({-1, -1, 0})
    }));
    par.idx_h = sys.size() - 1;
  }
};
#endif
