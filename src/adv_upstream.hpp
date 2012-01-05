/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 *  @section DESCRIPTION
 *    C++ implementation of the upstream/upwind/donor-cell scheme 
 *    for the Arakawa-C grid
 */
#ifndef ADV_UPSTREAM_HPP
#  define ADV_UPSTREAM_HPP

#  include "adv.hpp"
#  include "grd_arakawa-c-lorenz.hpp"

template <typename real_t> 
class adv_upstream : public adv<real_t> 
{
  adv_hack_macro // workaround for virtual template methods
  public: const int stencil_extent() { return 3; }
  public: const int time_levels() { return 2; }

  private: grd_arakawa_c_lorenz<real_t> *grid;

  public: adv_upstream(grd_arakawa_c_lorenz<real_t> *grid) 
    : grid(grid)
  { }

  public: 
  template <class idx>
  void op(
    arr<real_t>* psi[], 
    arr<real_t>* [], 
    arr<real_t>* tmp_v[], 
    const Range &i, const Range &j, const Range &k, 
    const int n, const int step,
    const arr<real_t> * const Cx, 
    const arr<real_t> * const, 
    const arr<real_t> * const 
  )
  {
    assert(step == 1);

    // eq. (2) in Smolarkiewicz & Margolin 1998

    /// \f$ F(\psi_l, \psi_r, U) = 0.5 \cdot (U + |U|) \cdot \psi_l + 0.5 \cdot (U - |U|) \cdot \psi_r \f$ \n
    /// eq. (3 a-d) in Smolarkiewicz & Margolin 1998 (J. Comp. Phys., 140, 459-480) \n
#  define mpdata_F(p1, p2, U) (real_t(.5) * (U + abs(U)) * p1 + real_t(.5) * (U - abs(U)) * p2)
    /// \f$ \psi^{n+1}_i = \psi^{n}_i - \left[ F( \psi^{n}_{i}, \psi^{n}_{i+1}, U_{i+1/2} ) - F( \psi^{n}_{i-1}, \psi^{n}_{i}, U_{i-1/2} ) \right] \f$ \n
    /// eq. (2) in Smolarkiewicz & Margolin 1998 (J. Comp. Phys., 140, 459-480) \n

    (*psi[n+1])(idx(i,j,k)) -= (
      mpdata_F((*psi[n])(idx(i  ,j,k)), (*psi[n])(idx(i+1,j,k)), (*Cx)(idx(i + grid->p_half,j,k))) - 
      mpdata_F((*psi[n])(idx(i-1,j,k)), (*psi[n])(idx(i  ,j,k)), (*Cx)(idx(i - grid->m_half,j,k)))
    );
#  undef mpdata_F
  }
};
#endif
