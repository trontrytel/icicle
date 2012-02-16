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
  private: vector<struct eqs<real_t>::gte> sys;
  public: vector<struct eqs<real_t>::gte> & system() { return sys; }

  private: quantity<si::length, real_t> h_unit;
  private: quantity<velocity_times_length, real_t> q_unit;
  public: eqs_shallow_water(const grd<real_t> &grid)
  {
    if (grid.nz() != 1) error_macro("only 1D (X or Y) and 2D (XY) simullations supported")

    h_unit = 1 * si::metres; 
    q_unit = 1 * si::metres * si::metres / si::seconds; 

    //if (grid.nx() != 1)
    {
      struct eqs<real_t>::gte e = {
        "qx", "heigh-integrated specific momentum (x)", 
        this->quan2str(q_unit), 
        1, 0, 0
      };
      sys.push_back(e);
    }
    //if (grid.ny() != 1) 
    {
      struct eqs<real_t>::gte e = {
        "qy", "heigh-integrated specific momentum (y)", 
        this->quan2str(q_unit), 
        0, 1, 0 
      };
      sys.push_back(e);
    }
    {
      struct eqs<real_t>::gte e = {
        "h", "thickness of the fluid layer", 
        this->quan2str(h_unit), 
        -1, -1, 0
      };
      sys.push_back(e);
    }
  }
};
#endif
