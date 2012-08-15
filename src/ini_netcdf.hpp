/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date January 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#pragma once
#include "ini.hpp"

#include <memory>
using std::unique_ptr; // TODO: detail

template <typename real_t>
class ini_netcdf : public ini<real_t>
{
  public: ini_netcdf(const grd<real_t> &grid, const string filename)
#if !defined(USE_NETCDF)
  {
     error_macro("recompile icicle with -DUSE_NETCDF")
  }
#else
  ;
  private: unique_ptr<NcFile> f;
#endif
  private: size_t xs, ys, zs;
  private: string filename;
  private: const grd<real_t> *grid;

  public: virtual void populate_scalar_field(
    const string &varname,
    const mtx::idx &ijk,
    mtx::arr<real_t> &data
  ) const;
};
