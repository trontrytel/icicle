/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#pragma once
#include "out.hpp"
#include "stp.hpp"
#include "cfg/cfg_netcdf.hpp"

#include <memory>  
using std::unique_ptr;

template <typename real_t>
class out_netcdf : public out<real_t>
{
  // ctor
  public: out_netcdf(
    const string &file, 
    const stp<real_t> &setup, 
    int ver, 
    const string &cmdline
  )
#if !defined(USE_NETCDF)
  { error_macro("recompile icicle with -DUSE_NETCDF (cmake -DNETCDF=ON)") }
#else
  ;

  // pimpl
  private: struct detail;
  private: unique_ptr<detail> pimpl;
#endif

  // dtor (e.g. writing execution time into the netcdf)
  public: ~out_netcdf();

  // implementing the public interface
  public: void record(
    const string &name, 
    const mtx::arr<real_t> &psi,
    const mtx::idx &ijk, 
    const unsigned long t // t is the number of the record!
  ); 
};
