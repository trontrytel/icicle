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

  private: int i_qx, i_qy;
  public: int idx_qx() { assert(i_qx >= 0); return i_qx; }
  public: int idx_qy() { assert(i_qy >= 0); return i_qy; }
  public: bool has_qx() { return i_qx != -1; }
  public: bool has_qy() { return i_qy != -1; }

  private: quantity<si::length, real_t> h_unit;
  private: quantity<si::momentum, real_t> q_unit;

  public: eqs_shallow_water(const grd<real_t> &grid)
    : i_qx(-1), i_qy(-1)
  {
    if (grid.nz() != 1) error_macro("only 1D (X) and 2D (XY) simullations supported")

    h_unit = 1 * si::metres;
    q_unit = 1 * si::newton * si::second;

    {
      struct eqs<real_t>::gte e = {"h", this->quan2str(h_unit)};
      sys.push_back(e);
    }
    if (grid.nx() != 1)
    {
      struct eqs<real_t>::gte e = {"qx", this->quan2str(q_unit)};
      sys.push_back(e);
      i_qx = sys.size() - 1;
    }
    if (grid.ny() != 1) 
    {
      struct eqs<real_t>::gte e = {"qy", this->quan2str(q_unit)};
      sys.push_back(e);
      i_qy = sys.size() - 1;
    }
cerr << "i_qx=" << i_qx << " i_qy=" << i_qy << endl;
  }

  public: vector<struct eqs<real_t>::gte> & system()
  {
    return sys;
  }
};
#endif
