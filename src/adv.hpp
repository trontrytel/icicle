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
    arr<real_t> *psi[], \
    arr<real_t> *tmp_s[], \
    arr<real_t> *tmp_v[], \
    const rng &i, const rng &j, const rng &k, \
    const int n, const int step, \
    const arr<real_t> * const Cx, const arr<real_t> * const Cy, const arr<real_t> * const Cz \
  ) { op<idx_ijk>(psi, tmp_s, tmp_v, i, j, k, n, step, Cx, Cy, Cz); } \
  public: void op_jki( \
    arr<real_t> *psi[], \
    arr<real_t> *tmp_s[], \
    arr<real_t> *tmp_v[], \
    const rng &i, const rng &j, const rng &k, \
    const int n, const int step, \
    const arr<real_t> * const Cx, const arr<real_t> * const Cy, const arr<real_t> * const Cz \
  ) { op<idx_jki>(psi, tmp_s, tmp_v, j, k, i, n, step, Cy, Cz, Cx); } \
  public: void op_kij( \
    arr<real_t> *psi[], \
    arr<real_t> *tmp_s[], \
    arr<real_t> *tmp_v[], \
    const rng &i, const rng &j, const rng &k, \
    const int n, const int step, \
    const arr<real_t> * const Cx, const arr<real_t> * const Cy, const arr<real_t> * const Cz \
  ) { op<idx_kij>(psi, tmp_s, tmp_v, k, i, j, n, step, Cz, Cx, Cy); }

template <typename real_t>
class adv : root
{
  public: virtual const int stencil_extent() = 0;
  public: virtual const int time_levels() = 0;
  public: virtual const int num_steps() { return 1; }
  public: virtual const int num_sclr_caches() { return 0; }
  public: virtual const int num_vctr_caches() { return 0; }

  private: virtual void op_ijk(
    arr<real_t> *psi[], 
    arr<real_t> *tmp_s[], 
    arr<real_t> *tmp_v[], 
    const rng &i, const rng &j, const rng &k, 
    const int n, const int step,
    const arr<real_t> * const Cx, const arr<real_t> * const Cy, const arr<real_t> * const Cz
  ) = 0;

  private: virtual void op_jki(
    arr<real_t> *psi[], 
    arr<real_t> *tmp_s[], 
    arr<real_t> *tmp_v[], 
    const rng &i, const rng &j, const rng &k, 
    const int n, const int step,
    const arr<real_t> * const Cx, const arr<real_t> * const Cy, const arr<real_t> * const Cz
  ) = 0; 

  private: virtual void op_kij(
    arr<real_t> *psi[], 
    arr<real_t> *tmp_s[], 
    arr<real_t> *tmp_v[], 
    const rng &i, const rng &j, const rng &k, 
    const int n, const int step,
    const arr<real_t> * const Cx, const arr<real_t> * const Cy, const arr<real_t> * const Cz
  ) = 0;

  public: virtual void op3D(
    arr<real_t> *psi[], 
    arr<real_t> *tmp_s[], 
    arr<real_t> *tmp_v[], 
    const rng &i, const rng &j, const rng &k, 
    const int n, const int s,
    const arr<real_t> * const Cx, const arr<real_t> * const Cy, const arr<real_t> * const Cz
  )
  {
    // op() uses the -= operator so the first assignment happens here
    *psi[n+1] = *psi[0];

    // we use the same code for each dimension switching the indices accordingly
    if (true)  
      op_ijk(psi, tmp_s, tmp_v, i, j, k, n, s, Cx, Cy, Cz); // X
    if (j.first() != j.last())
      op_jki(psi, tmp_s, tmp_v, i, j, k, n, s, Cx, Cy, Cz); // Y
    if (k.first() != k.last())
      op_kij(psi, tmp_s, tmp_v, i, j, k, n, s, Cx, Cy, Cz); // Z
  }
};
#endif
