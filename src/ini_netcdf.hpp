/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date January 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#pragma once

#include "ini.hpp"
#include "cfg/cfg_netcdf.hpp"

#include <memory>
using std::unique_ptr; 

template <typename real_t>
class ini_netcdf : public ini<real_t>
{
  // ctor
  public: ini_netcdf(const grd<real_t> &grid, const string filename)
#if !defined(USE_NETCDF)
  { error_macro("recompile icicle with -DUSE_NETCDF") }
#else
  ;
#endif

  // pimpl
  private: struct detail;
  private: unique_ptr<detail> pimpl;

  public: virtual void populate_scalar_field(
    const string &varname,
    const mtx::idx &ijk,
    mtx::arr<real_t> &data
  ) const;
};
