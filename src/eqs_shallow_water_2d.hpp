/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date January 2012
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef EQS_SHALLOW_WATER_2D_HPP
#  define EQS_SHALLOW_WATER_2D_HPP

#  include "cmn.hpp" // root class, error reporting
#  include "eqs.hpp"

template <typename real_t>
class eqs_shallow_water_2d : public eqs<real_t>
{
  private: vector<struct eqs<real_t>::gte> sys;

  private: quantity<si::length, real_t> h_unit;
  private: quantity<si::momentum, real_t> q_unit;

  public: eqs_shallow_water_2d()
  {
    h_unit = 1 * si::metres;
    q_unit = 1 * si::newton * si::second;

    {
      struct eqs<real_t>::gte e = {"h", this->quan2str(h_unit)};
      sys.push_back(e);
    }
    {
      struct eqs<real_t>::gte e = {"qx", this->quan2str(q_unit)};
      sys.push_back(e);
    }
    {
      struct eqs<real_t>::gte e = {"qy", this->quan2str(q_unit)};
      sys.push_back(e);
    }
  }

  public: vector<struct eqs<real_t>::gte> & system()
  {
    return sys;
  }
};
#endif
