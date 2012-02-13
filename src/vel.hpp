/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011 - February 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef VEL_HPP
#  define VEL_HPP

#  include "cmn.hpp" 
#  include "mtx.hpp" 
#  include "grd.hpp" 

template <typename real_t>
class vel : root
{
  public: virtual bool is_constant() = 0;

  public: virtual void populate_courant_fields(int n,
    mtx::arr<real_t> *Cx, 
    mtx::arr<real_t> *Cy, 
    mtx::arr<real_t> *Cz, 
    quantity<si::time, real_t> dt,
    mtx::arr<real_t> *Qx[],
    mtx::arr<real_t> *Qy[],
    mtx::arr<real_t> *Qz[]
  ) = 0;
};
#endif
