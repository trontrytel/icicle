/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef INI_HPP
#  define INI_HPP

#  include "grd.hpp"
#  include "mtx.hpp"

template <typename real_t>
class ini : root
{
  public: virtual void populate_scalar_field(
    const string &varname,
    const mtx::idx &ijk,
    mtx::arr<real_t> &data
  ) = 0;
};
#endif
