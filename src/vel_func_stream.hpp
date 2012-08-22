/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date April 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#pragma once
#include "vel_func.hpp" 

#include <map>
using std::map;

#include <string>
using std::string;

#include <vector>
using std::vector;

#include "cmn/cmn_netcdf.hpp" // TODO: move somehow into detail?
#include "cmn/cmn_error.hpp"

template <typename real_t>
class vel_func_stream : public vel_func<real_t>
{
  // stream-function works for 2D flows only
  public: virtual quantity<si::velocity, real_t> w(
    const quantity<si::length, real_t> &,  
    const quantity<si::length, real_t> &,  
    const quantity<si::length, real_t> &
  ) const
  {
    return real_t(0.) * si::metres / si::seconds; // TODO: assert(false);
  }

  private: map<quantity<si::length, real_t>, quantity<si::mass_density, real_t>> rhomap;
  
  protected: quantity<si::mass_density, real_t> rho(
    const quantity<si::length, real_t> &z
  ) const
  {
    assert(z >= 0 * si::metres);
    return rhomap.at(z);
  }

  public: vel_func_stream(const grd<real_t> &grid, const string &filename)
    : vel_func<real_t>(grid)
  { 
    // TODO: and what about USE_NETCDF=OFF !!!!
    try 
    {   
      NcFile nf(filename, NcFile::read);
      for (size_t j = 0; j < nf.getDim("Y").getSize(); ++j)
      {
        real_t tmp;
        nf.getVar("Y").getVar(vector<size_t>({j}), &tmp);
        quantity<si::length, real_t> z = tmp * si::metres;
        nf.getVar("rho").getVar(vector<size_t>({j}), &tmp);
        quantity<si::mass_density, real_t> rho = tmp * si::kilograms / si::cubic_metres;
        rhomap[z] = rho;
      }
    }
    catch (NcException& e) error_macro(e.what())
  }
};
