/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @brief definition of the adv_upstream class with an implementation of the upstream advection scheme
 */
#pragma once

#include "adv.hpp"
#include "grd.hpp"

/** @brief
 *  implementation of the upstream/upwind/donor-cell scheme 
 */
template <typename real_t> 
class adv_upstream : public adv<real_t> 
{
  public: const int stencil_extent() const { return 3; }
  public: const int time_levels() const { return 2; }

  // "Donor cell schemes are more accurate for larger Courant numbers"
  // Margolin & Smolarkiewicz 1998, SIAM J. Sci. Comput. (page 923)
  public: const real_t courant_max() const { return 1.; }
  public: const real_t courant_min() const { return .5; }

  public: class op3D : public adv<real_t>::op3D
  {
    // private nested class for storing indices
    private: 
    template <class idx>
    class indices
    {
      public: const idx i_j_k, iph_j_k, imh_j_k, ip1_j_k, im1_j_k;
      public: indices(
        const mtx::rng &i,
        const mtx::rng &j,
        const mtx::rng &k
      ) :
        i_j_k(  idx(i,               j, k)),
        iph_j_k(idx(i + static_rational<1,2>(), j, k)),
        imh_j_k(idx(i - static_rational<1,2>(), j, k)),
        ip1_j_k(idx(i + 1,           j, k)),
        im1_j_k(idx(i - 1,           j, k))
      { }
    };
 
    // member fields
    public: const indices<mtx::idx_ijk> indcs_x;
    public: const indices<mtx::idx_jki> indcs_y;
    public: const indices<mtx::idx_kij> indcs_z;

    // ctor
    public: op3D(const mtx::idx &ijk) :
      adv<real_t>::op3D(ijk),
      indcs_x(ijk.i, ijk.j, ijk.k),
      indcs_y(ijk.j, ijk.k, ijk.i), 
      indcs_z(ijk.k, ijk.i, ijk.j)
    { }

    // () operator
    public: virtual void operator()(
      mtx::arr<real_t> *psi[],
      const int n,
      const int s,
      const mtx::arr<real_t> * const Cx,
      const mtx::arr<real_t> * const Cy,
      const mtx::arr<real_t> * const Cz
    )
    { 
      *psi[n+1] = *psi[0]; // TODO: move from here to the solver?
      if (this->do_x()) op1D(psi, indcs_x, n, s, Cx, Cy, Cz);
      if (this->do_y()) op1D(psi, indcs_y, n, s, Cy, Cz, Cx);
      if (this->do_z()) op1D(psi, indcs_z, n, s, Cz, Cx, Cy); 
    }
    
    /// \f$ 
    ///   F(\psi_l, \psi_r, U) = 0.5 \cdot (U + |U|) \cdot \psi_l + 0.5 \cdot (U - |U|) \cdot \psi_r 
    /// \f$ 
    /// eq. (3 a-d) in Smolarkiewicz & Margolin 1998 (J. Comp. Phys., 140, 459-480)
    public: mtx_expr_3arg_macro(mpdata_F, p1, p2, U,
      adv<real_t>::pospart(U) * p1 + adv<real_t>::negpart(U) * p2
    )

    public: 
    template <class indices>
    void op1D(
      mtx::arr<real_t>* psi[], 
      const indices &idx,
      const int n,
      const int step,
      const mtx::arr<real_t> * const Cx,
      const mtx::arr<real_t> * const,
      const mtx::arr<real_t> * const
    )
    {
      assert(step == 1); 
// TODO: perpahs use static cast to get rid of the idx argument!

      /// \f$ 
      ///   \psi^{n+1}_i = \psi^{n}_i - \left[ 
      ///     F( \psi^{n}_{i}, \psi^{n}_{i+1}, U_{i+1/2} ) -
      ///     F( \psi^{n}_{i-1}, \psi^{n}_{i}, U_{i-1/2} ) 
      ///   \right] 
      /// \f$ 
      /// eq. (2) in Smolarkiewicz & Margolin 1998 (J. Comp. Phys., 140, 459-480) \n
      (*psi[n+1])(idx.i_j_k) -= (
        mpdata_F((*psi[n])(idx.i_j_k),   (*psi[n])(idx.ip1_j_k), (*Cx)(idx.iph_j_k)) - 
        mpdata_F((*psi[n])(idx.im1_j_k), (*psi[n])(idx.i_j_k),   (*Cx)(idx.imh_j_k))
      ); 

    };
  };

  public: typename adv<real_t>::op3D *factory(
    const mtx::idx &ijk,
    mtx::arr<real_t> **,
    mtx::arr<real_t> **,
    bool
  ) const;
};
