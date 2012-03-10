/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @brief definition of the adv_upstream class with an implementation of the upstream advection scheme
 */
#ifndef ADV_UPSTREAM_HPP
#  define ADV_UPSTREAM_HPP

#  include "adv.hpp"
#  include "grd_arakawa-c-lorenz.hpp"

/** @brief
 *  C++ implementation of the upstream/upwind/donor-cell scheme 
 *  for the Arakawa-C grid
 */
template <typename real_t> 
class adv_upstream : public adv<real_t> 
{
  public: const int stencil_extent() { return 3; }
  public: const int time_levels() { return 2; }

  // "Donor cell schemes are more accurate for larger Courant numbers"
  // Margolin & Smolarkiewicz 1998, SIAM J. Sci. Comput. (page 923)
  public: const real_t courant_max() { return 1.; }
  public: const real_t courant_min() { return .5; }

  private: grd_arakawa_c_lorenz<real_t> *grid;

  public: adv_upstream(grd_arakawa_c_lorenz<real_t> *grid) 
    : grid(grid)
  { }

  public: class op3D : public adv<real_t>::op3D
  {
    // private nested class for storing indices
    protected:
    template <class idx>
    class indices
    {
      public: idx i_j_k, iph_j_k, imh_j_k, ip1_j_k, im1_j_k;
      public: bool do_x, do_y, do_z;
      public: indices(
        const mtx::rng &i,
        const mtx::rng &j,
        const mtx::rng &k, 
        const grd_arakawa_c_lorenz<real_t> &grid
      ) :
        i_j_k(  idx(i,               j, k)),
        iph_j_k(idx(i + grid.p_half, j, k)),
        imh_j_k(idx(i - grid.m_half, j, k)),
        ip1_j_k(idx(i + 1,           j, k)),
        im1_j_k(idx(i - 1,           j, k)),
        do_x(i.first() != i.last()),
        do_y(j.first() != j.last()),
        do_z(k.first() != k.last())
      { }
    };

    // private members 
    private: indices<mtx::idx_ijk> idxx;
    private: indices<mtx::idx_jki> idxy;
    private: indices<mtx::idx_kij> idxz;

    // ctor
    public: op3D(const mtx::idx &ijk, const grd_arakawa_c_lorenz<real_t> &grid) :
      idxx(ijk.i, ijk.j, ijk.k, grid), 
      idxy(ijk.j, ijk.k, ijk.i, grid), 
      idxz(ijk.k, ijk.i, ijk.j, grid)
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
      if (idxx.do_x) op1D(psi, idxx, n, s, Cx, Cy, Cz);
      if (idxx.do_y) op1D(psi, idxy, n, s, Cy, Cz, Cx);
      if (idxx.do_z) op1D(psi, idxz, n, s, Cz, Cx, Cy); 
    }
    
    // protected methods
    // TODO: document, include as an option / autotest in CMake
#  ifdef MPDATA_NEGPOSPART_ABS
    protected: mtx_expr_1arg_macro(mpdata_pospart, x,
      real_t(.5) * (x + abs(x))
    )
    protected: mtx_expr_1arg_macro(mpdata_negpart, x,
      real_t(.5) * (x - abs(x))
    )
#  else
    protected: mtx_expr_1arg_macro(mpdata_pospart, x,
      max(real_t(0), x)
    )
    protected: mtx_expr_1arg_macro(mpdata_negpart, x,
      min(real_t(0), x)
    )
#  endif

    /// \f$ 
    ///   F(\psi_l, \psi_r, U) = 0.5 \cdot (U + |U|) \cdot \psi_l + 0.5 \cdot (U - |U|) \cdot \psi_r 
    /// \f$ 
    /// eq. (3 a-d) in Smolarkiewicz & Margolin 1998 (J. Comp. Phys., 140, 459-480)
    protected: mtx_expr_3arg_macro(mpdata_F, p1, p2, U,
      this->mpdata_pospart(U) * p1 + this->mpdata_negpart(U) * p2
    )

    protected: 
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
    }
  };

  public: virtual op3D *factory(
    const mtx::idx &ijk,
    mtx::arr<real_t> **,
    mtx::arr<real_t> **,
    bool
  )
  {
    return new op3D(ijk, *grid);
  }
};
#endif
