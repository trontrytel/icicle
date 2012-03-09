/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @brief definition of and helper macros for \ref adv - the base class for advection operators
 */
#ifndef ADV_HPP
#  define ADV_HPP

#  include "cmn.hpp" // root class, error reporting, ...
#  include "grd.hpp" // m_half, p_half, ...

/// @brief a base class for advection operators
template <typename real_t>
class adv : root
{
  public: virtual const int stencil_extent() = 0;
  public: virtual const int time_levels() = 0;
  public: virtual const int num_steps() { return 1; }
  public: virtual const int num_sclr_caches() { return 0; }
  public: virtual const int num_vctr_caches() { return 0; }

  public: virtual const real_t courant_max() = 0; 
  public: virtual const real_t courant_min() = 0;

  public: virtual void op3D(
    mtx::arr<real_t> *psi[], 
    mtx::arr<real_t> *tmp_s[], 
    mtx::arr<real_t> *tmp_v[], 
    const mtx::idx &ijk,
    const int n, const int s,
    const mtx::arr<real_t> * const Cx, const mtx::arr<real_t> * const Cy, const mtx::arr<real_t> * const Cz
  ) = 0;
};
#endif
