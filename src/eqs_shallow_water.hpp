/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date January 2012 - February 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef EQS_SHALLOW_WATER_HPP
#  define EQS_SHALLOW_WATER_HPP

#  include "cmn.hpp"
#  include "eqs.hpp"
#  include "grd.hpp"

template <typename real_t>
class eqs_shallow_water : public eqs<real_t> 
{
  private: ptr_vector<struct eqs<real_t>::gte> sys;
  public: ptr_vector<struct eqs<real_t>::gte> &system() { return sys; }

  private: class forcings : public rhs<real_t>
  {
    private: quantity<si::acceleration, real_t> g;
    private: quantity<si::length, real_t> dx;
    private: quantity<si::length, real_t> h_unit;
    private: quantity<velocity_times_length, real_t> q_unit;
    public: forcings(
      quantity<si::acceleration, real_t> g, 
      quantity<si::length, real_t> dx
    ) 
      : g(g), dx(dx)
    {
      h_unit = 1 * si::metres;  // TODO: repeated below
      q_unit = 1 * si::metres * si::metres / si::seconds;  // TODO: repeated below
    }

    // TODO: Coriolis
    // TODO: separate x and y forcings
    // TODO: 2 -> idx_h ?
    // TODO: stencil-extent-like method?
    public: void operator()(mtx::arr<real_t> &R, mtx::arr<real_t> **psi, mtx::idx &ijk) 
    { 
      R(ijk) -= real_t(g / (real_t(2) * dx) * h_unit * h_unit * si::seconds / q_unit) 
        * ((*psi[2])(ijk)) * 
        (
          ((*psi[2])(ijk.i + 1, ijk.j, ijk.k) /* TODO: + H0 */)
          - 
          ((*psi[2])(ijk.i - 1, ijk.j, ijk.k) /* TODO: + H0 */)
        );
    };
  };

  private: quantity<si::length, real_t> h_unit;
  private: quantity<velocity_times_length, real_t> q_unit;
  public: eqs_shallow_water(const grd<real_t> &grid)
  {
    if (grid.nz() != 1) error_macro("only 1D (X or Y) and 2D (XY) simullations supported")

    h_unit = 1 * si::metres; 
    q_unit = 1 * si::metres * si::metres / si::seconds; 

    //if (grid.nx() != 1) // TODO?
    {
      sys.push_back(new struct eqs<real_t>::gte({
        "qx", "heigh-integrated specific momentum (x)", 
        this->quan2str(q_unit), 
	vector<int>({1, 0, 0})
      }));
      sys.back().source_terms.push_back(new forcings(real_t(10.) * si::metres_per_second_squared, grid.dx())); // TODO ! 10
    }
    //if (grid.ny() != 1) // TODO?
    {
      sys.push_back(new struct eqs<real_t>::gte({
        "qy", "heigh-integrated specific momentum (y)", 
        this->quan2str(q_unit), 
        vector<int>({0, 1, 0})
      }));
      //sys.back().source_terms.push_back(new ...); // TODO!
    }
    {
      sys.push_back(new struct eqs<real_t>::gte({
        "h", "thickness of the fluid layer", 
        this->quan2str(h_unit), 
        vector<int>({-1, -1, 0})
      }));
    }
  }
};
#endif
