/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef ADV_HPP
#  define ADV_HPP

#  include "config.hpp" // USE_* defines
#  include "common.hpp" // root class, error reporting, ...
#  include "grd.hpp" // m_half, p_half, ...
#  include "idx.hpp"

template <typename real_t>
class adv : root
{
  public: virtual const int stencil_extent() = 0;
  public: virtual const int time_levels() = 0;
  public: virtual const int num_steps() { return 1; }
  public: virtual const int num_sclr_caches() { return 0; }
  public: virtual const int num_vctr_caches() { return 0; }

  private: virtual void op_ijk(
    Array<real_t, 3> *psi[], 
    Array<real_t, 3> *tmp_s[], 
    Array<real_t, 3> *tmp_v[], 
    const Range &i, const Range &j, const Range &k, 
    const int n, const int step,
    const Array<real_t, 3> &Cx, const Array<real_t, 3> &Cy, const Array<real_t, 3> &Cz
  ) = 0;

  private: virtual void op_jki(
    Array<real_t, 3> *psi[], 
    Array<real_t, 3> *tmp_s[], 
    Array<real_t, 3> *tmp_v[], 
    const Range &i, const Range &j, const Range &k, 
    const int n, const int step,
    const Array<real_t, 3> &Cx, const Array<real_t, 3> &Cy, const Array<real_t, 3> &Cz
  ) = 0; 

  private: virtual void op_kij(
    Array<real_t, 3> *psi[], 
    Array<real_t, 3> *tmp_s[], 
    Array<real_t, 3> *tmp_v[], 
    const Range &i, const Range &j, const Range &k, 
    const int n, const int step,
    const Array<real_t, 3> &Cx, const Array<real_t, 3> &Cy, const Array<real_t, 3> &Cz
  ) = 0;

  public: virtual void op3D(
    Array<real_t, 3> *psi[], 
    Array<real_t, 3> *tmp_s[], 
    Array<real_t, 3> *tmp_v[], 
    const Range &i, 
    const Range &j, 
    const Range &k, 
    const int n, const int s,
    Array<real_t, 3> &Cx,
    Array<real_t, 3> &Cy,
    Array<real_t, 3> &Cz
  )
  {
    // op() uses the -= operator so the first assignment happens here
    *psi[n+1] = *psi[0];

    if (true)  
      op_ijk(psi, tmp_s, tmp_v, i, j, k, n, s, Cx, Cy, Cz); // X
    if (j.first() != j.last())
      op_jki(psi, tmp_s, tmp_v, i, j, k, n, s, Cx, Cy, Cz); // Y
    if (k.first() != k.last())
      op_kij(psi, tmp_s, tmp_v, i, j, k, n, s, Cx, Cy, Cz); // Z
  }
};
#endif
