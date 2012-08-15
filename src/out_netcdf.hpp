/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#pragma once
#ifdef USE_NETCDF
#  include "out.hpp"
#  include "inf.hpp"
#  include "eqs.hpp"
#  include "stp.hpp"
#  include "grd.hpp"

#  include <memory>    // TODO: move into detail
using std::unique_ptr; // 

template <typename real_t>
class out_netcdf : public out<real_t>
{
  // TODO: move into detail
  private: unique_ptr<NcFile> f; 
  private: map<string, NcVar> vars;
  private: inf info;

  public: out_netcdf(
    const string &file, 
    const stp<real_t> &setup, 
    int ver, 
    const string &cmdline
  );

  public: ~out_netcdf();

  public: void record(
    const string &name, 
    const mtx::arr<real_t> &psi,
    const mtx::idx &ijk, 
    const unsigned long t // t is the number of the record!
  ); 
};
#endif 
