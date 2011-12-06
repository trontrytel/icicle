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
  public: virtual const int num_steps() = 0;

  public: virtual void op_ijk(
    Array<real_t, 3> *psi[], 
    const Range &i, const Range &j, const Range &k, 
    const int n, const int step,
    const Array<real_t, 3> &Cx, const Array<real_t, 3> &Cy, const Array<real_t, 3> &Cz
  ) = 0;

  public: virtual void op_jki(
    Array<real_t, 3> *psi[], 
    const Range &i, const Range &j, const Range &k, 
    const int n, const int step,
    const Array<real_t, 3> &Cx, const Array<real_t, 3> &Cy, const Array<real_t, 3> &Cz
  ) = 0; 

  public: virtual void op_kij(
    Array<real_t, 3> *psi[], 
    const Range &i, const Range &j, const Range &k, 
    const int n, const int step,
    const Array<real_t, 3> &Cx, const Array<real_t, 3> &Cy, const Array<real_t, 3> &Cz
  ) = 0;

/*
  public: virtual void op3D(
    Array<real_t, 3> *psi_ijk[], 
    Array<real_t, 3> *psi_jki[], 
    Array<real_t, 3> *psi_kij[], 
    const Range &i, 
    const Range &j, 
    const Range &k, 
    const int n, const int step,
    arr<real_t> &Cx,
    arr<real_t> &Cy,
    arr<real_t> &Cz
  )
  {
    // op() uses the -= operator so the first assignment happens here
    *psi_ijk[n+1] = *psi_ijk[0];

    if (true) // in extreme cases paralellisaion may make i->first() = i->last()
      op(0, psi_ijk, i, j, k, n, step, Cx.ijk(), Cy.ijk(), Cz.ijk()); // X
    if (j.first() != j.last())
      op(1, psi_jki, j, k, i, n, step, Cy.jki(), Cz.jki(), Cx.jki()); // Y
    if (k.first() != k.last())
      op(2, psi_kij, k, i, j, n, step, Cz.kij(), Cx.kij(), Cy.kij()); // Z
  }
*/
};
#endif
