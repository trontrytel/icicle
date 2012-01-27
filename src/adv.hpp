/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef ADV_HPP
#  define ADV_HPP

#  include "cmn.hpp" // root class, error reporting, ...
#  include "grd.hpp" // m_half, p_half, ...

// C++ forbins virtual template methods so we do a preprocessor workaround:
// the macro below has to be included in every class inheriting from adv
#  define adv_hack_macro \
  public: void op_ijk( \
    mtx::arr<real_t> *psi[], \
    mtx::arr<real_t> *tmp_s[], \
    mtx::arr<real_t> *tmp_v[], \
    const mtx::rng &i, const mtx::rng &j, const mtx::rng &k, \
    const int n, const int step, \
    const mtx::arr<real_t> * const Cx, const mtx::arr<real_t> * const Cy, const mtx::arr<real_t> * const Cz \
  ) { op<mtx::idx_ijk>(psi, tmp_s, tmp_v, i, j, k, n, step, Cx, Cy, Cz); } \
  public: void op_jki( \
    mtx::arr<real_t> *psi[], \
    mtx::arr<real_t> *tmp_s[], \
    mtx::arr<real_t> *tmp_v[], \
    const mtx::rng &i, const mtx::rng &j, const mtx::rng &k, \
    const int n, const int step, \
    const mtx::arr<real_t> * const Cx, const mtx::arr<real_t> * const Cy, const mtx::arr<real_t> * const Cz \
  ) { op<mtx::idx_jki>(psi, tmp_s, tmp_v, j, k, i, n, step, Cy, Cz, Cx); } \
  public: void op_kij( \
    mtx::arr<real_t> *psi[], \
    mtx::arr<real_t> *tmp_s[], \
    mtx::arr<real_t> *tmp_v[], \
    const mtx::rng &i, const mtx::rng &j, const mtx::rng &k, \
    const int n, const int step, \
    const mtx::arr<real_t> * const Cx, const mtx::arr<real_t> * const Cy, const mtx::arr<real_t> * const Cz \
  ) { op<mtx::idx_kij>(psi, tmp_s, tmp_v, k, i, j, n, step, Cz, Cx, Cy); }

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

  private: virtual void op_ijk(
    mtx::arr<real_t> *psi[], 
    mtx::arr<real_t> *tmp_s[], 
    mtx::arr<real_t> *tmp_v[], 
    const mtx::rng &i, const mtx::rng &j, const mtx::rng &k, 
    const int n, const int step,
    const mtx::arr<real_t> * const Cx, const mtx::arr<real_t> * const Cy, const mtx::arr<real_t> * const Cz
  ) = 0;

  private: virtual void op_jki(
    mtx::arr<real_t> *psi[], 
    mtx::arr<real_t> *tmp_s[], 
    mtx::arr<real_t> *tmp_v[], 
    const mtx::rng &i, const mtx::rng &j, const mtx::rng &k, 
    const int n, const int step,
    const mtx::arr<real_t> * const Cx, const mtx::arr<real_t> * const Cy, const mtx::arr<real_t> * const Cz
  ) = 0; 

  private: virtual void op_kij(
    mtx::arr<real_t> *psi[], 
    mtx::arr<real_t> *tmp_s[], 
    mtx::arr<real_t> *tmp_v[], 
    const mtx::rng &i, const mtx::rng &j, const mtx::rng &k, 
    const int n, const int step,
    const mtx::arr<real_t> * const Cx, const mtx::arr<real_t> * const Cy, const mtx::arr<real_t> * const Cz
  ) = 0;

  public: virtual void op3D(
    mtx::arr<real_t> *psi[], 
    mtx::arr<real_t> *tmp_s[], 
    mtx::arr<real_t> *tmp_v[], 
    const mtx::idx &ijk,
    const int n, const int s,
    const mtx::arr<real_t> * const Cx, const mtx::arr<real_t> * const Cy, const mtx::arr<real_t> * const Cz
  )
  {
    // op() uses the -= operator so the first assignment happens here
    *psi[n+1] = *psi[0];

    // we use the same code for each dimension switching the indices accordingly
    if (ijk.i_spans) op_ijk(psi, tmp_s, tmp_v, ijk.i, ijk.j, ijk.k, n, s, Cx, Cy, Cz); // X
    if (ijk.j_spans) op_jki(psi, tmp_s, tmp_v, ijk.i, ijk.j, ijk.k, n, s, Cx, Cy, Cz); // Y
    if (ijk.k_spans) op_kij(psi, tmp_s, tmp_v, ijk.i, ijk.j, ijk.k, n, s, Cx, Cy, Cz); // Z
  }
};
#endif
