/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011 - January 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @section DESCRIPTION
 *    C++ implementation of the second- and third-order MPDATA scheme 
 *    for solenoidal flows of scalar fields on a uniformly spaced 
 *    1D, 2D and 3D Arakawa-C grid 
 */
#ifndef ADV_MPDATA_HPP
#  define ADV_MPDATA_HPP

#  include "adv_upstream.hpp"
#  include "grd_arakawa-c-lorenz.hpp"

template <typename real_t> 
class adv_mpdata : public adv_upstream<real_t> 
{
  public: const int stencil_extent() 
  {
    int halo = 1;
    if (cross_terms && iord > 2) halo += (iord - 2);
    if (third_order && halo < 2) halo = 2;
    return 1 + 2 * halo;
  }
  public: const int num_steps() { return iord; }
  public: const int num_vctr_caches() 
  { 
    if (iord == 1) return 0;
    if (iord == 2) return 1; // storing anti-diff vel. instead of calculating it twice: for u_{i+1/2} and for u_{i-1/2}
    if (iord == 3) return 4; // storing all components of anti-diff vel. for the next iteration + ditto
    return 6; // for higher iords we have to save all components repeatedly
  }

  private: int iord;
  private: bool cross_terms, third_order;
  private: grd_arakawa_c_lorenz<real_t> *grid;

  public: adv_mpdata(grd_arakawa_c_lorenz<real_t> *grid, int iord, bool cross_terms, bool third_order) // TODO: enums?
    : iord(iord), cross_terms(cross_terms), third_order(third_order), grid(grid), adv_upstream<real_t>(grid)
  {
    if (iord <= 0) error_macro("iord (the number of iterations) must be > 0")
  }

  // TODO: enclose all arguments in parenthesis, i.e. U -> (U)
  // using preprocessor macros as it's tricky make methods return parts of Blitz expressions 
#    ifdef MPDATA_FRAC_EPSILON
#      define mpdata_frac(num, den) ((num) / (den + mtx::eps<real_t>()))
#    else
#      define mpdata_frac(num, den) (where(den > real_t(0), (num) / (den), real_t(0)))
#    endif

  // macros for 2nd order terms:
#    define mpdata_A(pr, pl) mpdata_frac(pr - pl, pr + pl) 
#    define mpdata_B(pru, plu, prd, pld) (real_t(.5) * mpdata_frac(pru + plu - prd - pld, pru + plu + prd + pld))
#    define mpdata_V(Vru, Vlu, Vrd, Vld) (real_t(.25) * (Vru + Vlu + Vrd + Vld))
#    define mpdata_W(Wru, Wlu, Wrd, Wld) mpdata_V(Wru, Wlu, Wrd, Wld)
#    define mpdata_CA(pr, pl, U) ((abs(U) - pow(U,2)) * mpdata_A(pr, pl))
#    define mpdata_CB(pru, plu, prd, pld, U, V) (U * V * mpdata_B(pru, plu, prd, pld)) 
  // macros for 3rd order terms:

  /// first term from eq. (36) from Smolarkiewicz & Margolin 1998 (with G=1)
  /// \f$ 
  ///   \frac{(\delta x)^2}{6} \left( 3U|U| - 3U^3 - U \right) 
  ///   \frac{1}{\psi} \frac{\partial^2 \psi}{\partial x^2} \approx
  ///   (\ldots) \cdot \frac{1}{3} \cdot 
  ///   \frac{\psi_{i+1} - \psi_{i+1} - \psi{i} + \psi{i-1}}{\psi_{i+1} + \psi_{i+1} + \psi{i} + \psi{i-1}}
  /// \f$ \n 
#    define mpdata_3rd_xx(pp2, pp1, pp0, pm1, U) (\
       (real_t(3) * U * abs(U) - real_t(2) * U * pow(U,3) - U) / real_t(3) \
       * mpdata_frac(pp2 - pp1 - pp0 + pm1, pp2 + pp1 + pp0 + pm1) \
     )
  /// second term from eq. (36) from Smolarkiewicz & Margolin 1998 (with G=1)
  /// \f$
  ///   \frac{\delta x \delta y}{2} V \left( |U| - 2U^2 \right)
  ///   \frac{1}{\psi} \frac{\partial^2 \psi}{\partial x \partial y} \approx
  ///   V (|U| - 2U^2) \frac{
  ///     \psi_{i+1,j+1} - \psi_{i,j+1} - \psi_{i+1,j-1} + \psi_{i,j-1}
  ///   }{
  ///     \psi_{i+1,j+1} + \psi_{i,j+1} + \psi_{i+1,j-1} + \psi_{i,j-1}
  ///   }
  /// \f$
#    define mpdata_3rd_xy(pip0jp1, pip0jm1, pip1jp1, pip1jm1, U, V) (\
       V * (abs(U) - real_t(2) * pow(U,2)) \
       * mpdata_frac( \
         pip1jp1 - pip0jp1 - pip1jm1 + pip0jm1, \
         pip1jp1 + pip0jp1 + pip1jm1 + pip0jm1 \
       ) \
     )
  /// third term from eq. (36) from Smolarkiewicz & Margolin 1998 (with G=1)
#    define mpdata_3rd_xz(pip0kp1, pip0km1, pip1kp1, pip1km1, U, W) \
       mpdata_3rd_xy(pip0kp1, pip0km1, pip1kp1, pip1km1, U, W)
  /// fourth term from eq. (36) from Smolarkiewicz & Margolin 1998 (with G=1)
  /// \f$
  ///   \frac{2}{3} \delta x \delta z UVW \frac{1}{\psi} \frac{\partial^2 \psi}{\partial y \partial z} \approx
  ///   \frac{2}{3} UVW \frac{
  ///     + \psi_{i,j+1,k+1} - \psi_{i,j-1,k+1} + \psi_{i+1,j+1,k+1} - \psi_{i+1,j-1,k+1} 
  ///     - \psi_{i,j+1,k-1} + \psi_{i,j-1,k-1} - \psi_{i+1,j+1,k-1} + \psi_{i+1,j-1,k-1}
  ///   }{
  ///     + \psi_{i,j+1,k+1} + \psi_{i,j-1,k+1} + \psi_{i+1,j+1,k+1} + \psi_{i+1,j-1,k+1} 
  ///     + \psi_{i,j+1,k-1} + \psi_{i,j-1,k-1} + \psi_{i+1,j+1,k-1} + \psi_{i+1,j-1,k-1}
  ///   }
  /// \f$
#    define mpdata_3rd_yz( \
       pip0jp1kp1, pip0jm1kp1, pip1jp1kp1, pip1jm1kp1, \
       pip0jp1km1, pip0jm1km1, pip1jp1km1, pip1jm1km1, \
       U, V, W \
     ) ( \
       real_t(2./3.) * U * V * W * mpdata_frac( \
         pip0jp1kp1 - pip0jm1kp1 + pip1jp1kp1 - pip1jm1kp1 - \
         pip0jp1km1 + pip0jm1km1 - pip1jp1km1 + pip1jm1km1, \
         pip0jp1kp1 + pip0jm1kp1 + pip1jp1kp1 + pip1jm1kp1 + \
         pip0jp1km1 + pip0jm1km1 + pip1jp1km1 + pip1jm1km1 \
       ) \
     )

  protected:
  template <class idx>
  void mpdata_U(
    const mtx::arr<real_t> * C_adf,
    const mtx::arr<real_t> * const psi[], const int n, const int step,
    const mtx::rng &i, const mtx::rng &j, const mtx::rng &k,
    const mtx::arr<real_t> &Cx, const mtx::arr<real_t> &Cy, const mtx::arr<real_t> &Cz
  )
  {
    /// multidimensional antidiffusive velocity: \n
    /// \f$ \tilde{U}^{I}_{i+1/2}=\left[ |U^{I}_{i+1/2}| \Delta x^{I} - \Delta t (u^{I}_{i+1/2})^2 \right] \cdot
    ///  \frac{\psi^{*}_{i+1}-\psi^{*}_{i}}{(\psi^{*}_{i+1}+\psi^{*}_{i}) \Delta x^{I}} -  \f$ \n
    /// \f$ - \sum\limits_{J=1,J \ne I} 0.5 \Delta t U^{I}_{i+1/2} \bar{U}^{J}_{i+1/2} \cdot
    /// \frac{\psi^{*}_{i+1,j+1}+\psi^{*}_{i,j+1}-\psi^{*}_{i+1,j-1}-\psi^{*}_{i,j-1}}
    /// {\psi^{*}_{i+1,j+1}+\psi^{*}_{i,j+1}+\psi^{*}_{i+1,j-1}+\psi^{*}_{i,j-1}} \f$ \n
    /// eq. (13-14) in Smolarkiewicz 1984 (J. Comp. Phys.,54,352-362) \n

    int iord_halo_yz = (iord > 2 && cross_terms) ? iord - step : 0; 
    int iord_halo_x  = (iord > 3 && cross_terms && iord != step) ? iord - step - 1 : 0;
    mtx::rng // modified indices
      im(i.first() - 1 - iord_halo_x, i.last() + iord_halo_x), // instead of computing u_{i+1/2} and u_{i-1/2} for all i we compute u_{i+1/2} for im=(i-1, ... i)
      jm(j.first() - iord_halo_yz, j.last() + iord_halo_yz), // 
      km(k.first() - iord_halo_yz, k.last() + iord_halo_yz); // 
    mtx::rng // forward-in-space perspective
      ir = im + 1,            // right
      ic = im + grid->p_half, // center
      il = im;                // left
    idx // output indices
      adfidx = idx(
        mtx::rng(i.first() - grid->m_half - iord_halo_x, i.last() + grid->p_half + iord_halo_x), 
        jm, 
        km
    );

    assert(finite(sum((*psi[n])(idx(ir, jm, km)))));
    assert(finite(sum((*psi[n])(idx(il, jm, km)))));
    assert(finite(sum(Cx(idx(ic, jm, km)))));

    (*C_adf)(adfidx) = (
      mpdata_CA( 
        (*psi[n])(idx(ir, jm, km)), (*psi[n])(idx(il, jm, km)), /* pl, pr */ 
        Cx(idx(ic, jm, km)) 
      )
    );  
    if (cross_terms) 
    {
      if (j.first() != j.last()) 
      {
        assert(finite(sum((*psi[n])(idx(ir, jm+1, km)))));
        assert(finite(sum((*psi[n])(idx(il, jm+1, km)))));
        assert(finite(sum((*psi[n])(idx(ir, jm-1, km)))));
        assert(finite(sum((*psi[n])(idx(il, jm-1, km)))));
        assert(finite(sum(Cy(idx(ir, jm + grid->p_half, km)))));
        assert(finite(sum(Cy(idx(il, jm + grid->p_half, km)))));
        assert(finite(sum(Cy(idx(ir, jm - grid->m_half, km)))));
        assert(finite(sum(Cy(idx(il, jm - grid->m_half, km)))));

        (*C_adf)(adfidx) -= (
          mpdata_CB( 
            (*psi[n])(idx(ir, jm+1, km)), (*psi[n])(idx(il, jm+1, km)), /* pru, plu */ 
            (*psi[n])(idx(ir, jm-1, km)), (*psi[n])(idx(il, jm-1, km)), /* prd, pld */ 
            Cx(idx(ic, jm, km)), mpdata_V( 
              Cy(idx(ir, jm + grid->p_half, km)), Cy(idx(il, jm + grid->p_half, km)), /* Vru, Vlu */ 
              Cy(idx(ir, jm - grid->m_half, km)), Cy(idx(il, jm - grid->m_half, km))  /* Vrd, Vld */ 
            ) 
          )
        );
      }
      if (k.first() != k.last()) 
      {
        assert(finite(sum((*psi[n])(idx(ir, jm, km+1)))));
        assert(finite(sum((*psi[n])(idx(il, jm, km+1)))));
        assert(finite(sum((*psi[n])(idx(ir, jm, km-1)))));
        assert(finite(sum((*psi[n])(idx(il, jm, km-1)))));
        assert(finite(sum(Cz(idx(ir, jm, km + grid->p_half)))));
        assert(finite(sum(Cz(idx(il, jm, km + grid->p_half)))));
        assert(finite(sum(Cz(idx(ir, jm, km - grid->m_half)))));
        assert(finite(sum(Cz(idx(il, jm, km - grid->m_half)))));

        (*C_adf)(adfidx) -= ( // otherwise Cz is uninitialised!
          mpdata_CB( 
            (*psi[n])(idx(ir, jm, km+1)), (*psi[n])(idx(il, jm, km+1)), /* pru, plu */ 
            (*psi[n])(idx(ir, jm, km-1)), (*psi[n])(idx(il, jm, km-1)), /* prd, pld */ 
            Cx(idx(ic, jm, km)), mpdata_W( 
              Cz(idx(ir, jm, km + grid->p_half)), Cz(idx(il, jm, km + grid->p_half)), /* Wru, Wlu */ 
              Cz(idx(ir, jm, km - grid->m_half)), Cz(idx(il, jm, km - grid->m_half))  /* Wrd, Wld */ 
            ) 
          )  
        );
      }
      if (third_order) 
      {
        if (j.first() != j.last()) 
        {
          (*C_adf)(adfidx) += ( // otherwise Cy is uninitialised
            mpdata_3rd_xy(
              (*psi[n])(idx(im  , jm+1, km)), // pip0jp1, 
              (*psi[n])(idx(im  , jm-1, km)), // pip0jm1, 
              (*psi[n])(idx(im+1, jm+1, km)), // pip1jp1, 
              (*psi[n])(idx(im+1, jm-1, km)), // pip1jm1, 
              Cx(idx(ic,jm,km)), // U, 
              mpdata_V( // V
                Cy(idx(ir, jm + grid->p_half, km)), Cy(idx(il, jm + grid->p_half, km)), /* Vru, Vlu */ 
                Cy(idx(ir, jm - grid->m_half, km)), Cy(idx(il, jm - grid->m_half, km))  /* Vrd, Vld */ 
              ) 
            )
          );
        }
        if (k.first() != k.last()) 
        {
          (*C_adf)(adfidx) += ( // otherwise Cz is uninitialised
            mpdata_3rd_xz(
              (*psi[n])(idx(im  , jm, km+1)), // pip0kp1, 
              (*psi[n])(idx(im  , jm, km-1)), // pip0km1, 
              (*psi[n])(idx(im+1, jm, km+1)), // pip1kp1, 
              (*psi[n])(idx(im+1, jm, km-1)), // pip1km1, 
              Cx(idx(ic,jm,km)), // U, 
              mpdata_W( // W
                Cz(idx(ir, jm, km + grid->p_half)), Cz(idx(il, jm, km + grid->p_half)), /* Wru, Wlu */ 
                Cz(idx(ir, jm, km - grid->m_half)), Cz(idx(il, jm, km - grid->m_half))  /* Wrd, Wld */ 
              ) 
            )
          ); 
        }
        if ((j.first() != j.last()) && (k.first() != k.last())) 
        {
          assert(finite(sum((*psi[n])(idx(im  , jm+1, km+1)))));
          assert(finite(sum((*psi[n])(idx(im  , jm-1, km+1)))));
          assert(finite(sum((*psi[n])(idx(im+1, jm+1, km+1)))));
          assert(finite(sum((*psi[n])(idx(im+1, jm-1, km+1)))));
          assert(finite(sum((*psi[n])(idx(im  , jm+1, km-1)))));
          assert(finite(sum((*psi[n])(idx(im  , jm-1, km-1)))));
          assert(finite(sum((*psi[n])(idx(im+1, jm+1, km-1)))));
          assert(finite(sum((*psi[n])(idx(im+1, jm-1, km-1)))));

          (*C_adf)(adfidx) -= ( // otherwise Cx & Cz are uninitialised
            mpdata_3rd_yz(
              (*psi[n])(idx(im  , jm+1, km+1)), // pip0jp1kp1, 
              (*psi[n])(idx(im  , jm-1, km+1)), // pip0jm1kp1, 
              (*psi[n])(idx(im+1, jm+1, km+1)), // pip1jp1kp1, 
              (*psi[n])(idx(im+1, jm-1, km+1)), // pip1jm1kp1, 
              (*psi[n])(idx(im  , jm+1, km-1)), // pip0jp1km1, 
              (*psi[n])(idx(im  , jm-1, km-1)), // pip0jm1km1, 
              (*psi[n])(idx(im+1, jm+1, km-1)), // pip1jp1km1, 
              (*psi[n])(idx(im+1, jm-1, km-1)), // pip1jm1km1, 
              Cx(idx(ic,jm,km)), // U
              mpdata_V( // V
                Cy(idx(ir, jm + grid->p_half, km)), Cy(idx(il, jm + grid->p_half, km)), /* Vru, Vlu */ 
                Cy(idx(ir, jm - grid->m_half, km)), Cy(idx(il, jm - grid->m_half, km))  /* Vrd, Vld */ 
              ), 
              mpdata_W( // W
                Cz(idx(ir, jm, km + grid->p_half)), Cz(idx(il, jm, km + grid->p_half)), /* Wru, Wlu */ 
                Cz(idx(ir, jm, km - grid->m_half)), Cz(idx(il, jm, km - grid->m_half))  /* Wrd, Wld */ 
              ) 
            )
          );
        }
      }
    }
    if (third_order)
    { 
      assert(finite(sum((*psi[n])(idx(im+2, jm, km)))));
      assert(finite(sum((*psi[n])(idx(im-1, jm, km)))));

      (*C_adf)(adfidx) +=(
        mpdata_3rd_xx(
          (*psi[n])(idx(im+2, jm, km)), 
          (*psi[n])(idx(im+1, jm, km)), 
          (*psi[n])(idx(im  , jm, km)), 
          (*psi[n])(idx(im-1, jm, km)), 
          Cx(idx(ic,jm,km))
        )
      );  
    }
  }

    /// 
    /// \f$ \psi^{n+1}_{i} = \psi^{*}_{i} -\sum\limits_{I} \left[ F^{I}( \psi^{*}_{i}, \psi^{*}_{i+1}, \tilde{U}^{I}_{i+1/2} )
    /// - F^{I}( \psi^{*}_{i-1}, \psi^{*}_{i}, \tilde{U}^{I}_{i-1/2} ) \right] \f$ \n
    /// where \f$ I \f$ denotes the sum over all dimensions and \f$ \tilde{U} \f$ is multidimensional antidiffusive velocity \n
    /// eq. (12) in Smolarkiewicz 1984 (J. Comp. Phys.,54,352-362) 
  public: void op3D(
    mtx::arr<real_t> *psi[], mtx::arr<real_t> *[], mtx::arr<real_t> *tmp_v[],
    const mtx::idx &ijk,
    const int n, const int step,
    const mtx::arr<real_t> * const Cx, const mtx::arr<real_t> * const Cy, const mtx::arr<real_t> * const Cz
  )
  {
    int 
      x_old = -1, y_old = -1, z_old = -1, // \__valid for iord < 3
      x_new =  0, y_new =  0, z_new =  0; // /  
    if (iord > 2)
    {
      x_old = 0, y_old = 1, z_old = 2;
      x_new = 0, y_new = 1, z_new = 2;
      if (iord == 3 && step == 3) x_new = y_new = z_new = 3; // saving memory for iord=3
      else
      {
        if (step % 2 == 0) 
          x_old += 3, y_old += 3, z_old += 3;
        else 
          x_new += 3, y_new += 3, z_new += 3;
      }
    }

    const mtx::arr<real_t> 
      * const Cx_unco = (step < 3 ? Cx : tmp_v[x_old]), 
      *       Cx_corr = (step < 2 ? Cx : tmp_v[x_new]),
      * const Cy_unco = (step < 3 ? Cy : tmp_v[y_old]), 
      *       Cy_corr = (step < 2 ? Cy : tmp_v[y_new]),
      * const Cz_unco = (step < 3 ? Cz : tmp_v[z_old]), 
      *       Cz_corr = (step < 2 ? Cz : tmp_v[z_new]);
     
    *psi[n+1] = *psi[0]; // TODO: at least this should be placed in adv... and the leapfrog & upstream in adv_dimsplit?
    if (ijk.i_spans)
    {
      if (step > 1) mpdata_U<mtx::idx_ijk>(Cx_corr, psi, n, step, ijk.i, ijk.j, ijk.k, *Cx_unco, *Cy_unco, *Cz_unco);
      adv_upstream<real_t>::template op<mtx::idx_ijk>(psi, NULL, NULL, ijk.i, ijk.j, ijk.k, n, 1, Cx_corr, NULL, NULL);
    }
    if (ijk.j_spans)
    {
      if (step > 1) mpdata_U<mtx::idx_jki>(Cy_corr, psi, n, step, ijk.j, ijk.k, ijk.i, *Cy_unco, *Cz_unco, *Cx_unco);
      adv_upstream<real_t>::template op<mtx::idx_jki>(psi, NULL, NULL, ijk.j, ijk.k, ijk.i, n, 1, Cy_corr, NULL, NULL);
    }
    if (ijk.k_spans)
    {
      if (step > 1) mpdata_U<mtx::idx_kij>(Cz_corr, psi, n, step, ijk.k, ijk.i, ijk.j, *Cz_unco, *Cx_unco, *Cy_unco);
      adv_upstream<real_t>::template op<mtx::idx_kij>(psi, NULL, NULL, ijk.k, ijk.i, ijk.j, n, 1, Cz_corr, NULL, NULL);
    }
  }
};
#endif
