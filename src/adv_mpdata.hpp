/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 *  @section DESCRIPTION
 *    C++ implementation of the MPDATA scheme for the Arakawa-C grid
 *    (for solenoidal flows on a uniformly spaced grid)
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
  public: virtual const int num_vctr_caches() 
  { 
    if (iord == 1) return 0;
    if (iord == 2) return 1; // storing anti-diff vel. instead of calculating it twice: for u_{i+1/2} and for u_{i-1/2}
    if (iord == 3) return 4; // storing all components of anti-diff vel. for the next iteration + ditto
    return 6; // for higher iords we have to save all components repeatedly
  }

  private: int iord;
  // TODO; cross_term = !1D
  private: bool cross_terms, third_order;
  private: grd_arakawa_c_lorenz<real_t> *grid;

  public: adv_mpdata(grd_arakawa_c_lorenz<real_t> *grid, int iord, bool cross_terms, bool third_order) // TODO: enums?
    : iord(iord), cross_terms(cross_terms), third_order(third_order), grid(grid), adv_upstream<real_t>(grid)
  {
    if (iord <= 0) error_macro("iord (the number of iterations) must be > 0")
  }

  // TODO: enclose all arguments in parenthesis, i.e. U -> (U)
  // using preprocessor macros as it's tricky make methods return parts of Blitz expressions 
#    define mpdata_frac(num, den) (where(den > real_t(0), (num) / (den), real_t(0)))
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
    const arr<real_t> * C_adf,
    arr<real_t> *psi[], const int n, const int step,
    const rng &i, const rng &j, const rng &k,
    const arr<real_t> &Cx, const arr<real_t> &Cy, const arr<real_t> &Cz
  )
  {
    /// multidimensional antidiffusive velocity: \n
    /// \f$ \tilde{U}^{I}_{i+1/2}=\left[ |U^{I}_{i+1/2}| \Delta x^{I} - \Delta t (u^{I}_{i+1/2})^2 \right] \cdot
    ///  \frac{\psi^{*}_{i+1}-\psi^{*}_{i}}{(\psi^{*}_{i+1}+\psi^{*}_{i}) \Delta x^{I}} -  \f$ \n
    /// \f$ - \sum\limits_{J=1,J \ne I} 0.5 \Delta t U^{I}_{i+1/2} \bar{U}^{J}_{i+1/2} \cdot
    /// \frac{\psi^{*}_{i+1,j+1}+\psi^{*}_{i,j+1}-\psi^{*}_{i+1,j-1}-\psi^{*}_{i,j-1}}
    /// {\psi^{*}_{i+1,j+1}+\psi^{*}_{i,j+1}+\psi^{*}_{i+1,j-1}+\psi^{*}_{i,j-1}} \f$ \n
    /// eq. (13-14) in Smolarkiewicz 1984 (J. Comp. Phys.,54,352-362) \n

    int iord_halo_yz = (iord > 2 && cross_terms) ? iord - step     : 0; 
    int iord_halo_x  = (iord > 3 && cross_terms) ? iord - step - 1 : 0;
    rng // modified indices
      im(i.first() - 1 - iord_halo_x, i.last() + iord_halo_x), // instead of computing u_{i+1/2} and u_{i-1/2} for all i we compute u_{i+1/2} for im=(i-1, ... i)
      jm(j.first() - iord_halo_yz, j.last() + iord_halo_yz), // 
      km(k.first() - iord_halo_yz, k.last() + iord_halo_yz); // 
    rng // forward-in-space perspective
      ir = im + 1,            // right
      ic = im + grid->p_half, // center
      il = im;                // left
    idx // output indices
      adfidx = idx(
        rng(i.first() - grid->m_half - iord_halo_x, i.last() + grid->p_half + iord_halo_x), 
        jm, 
        km
    );

    (*C_adf)(adfidx) = (
      mpdata_CA( 
        (*psi[n])(idx(ir, jm, km)), (*psi[n])(idx(il, jm, km)), /* pl, pr */ 
        Cx(idx(ic, jm, km)) 
      )
    );  
    if (cross_terms) 
    {
      (*C_adf)(adfidx) -= (
        mpdata_CB( 
          (*psi[n])(idx(ir, jm+1, km)), (*psi[n])(idx(il, jm+1, km)), /* pru, plu */ 
          (*psi[n])(idx(ir, jm-1, km)), (*psi[n])(idx(il, jm-1, km)), /* prd, pld */ 
          Cx(idx(ic, jm, km)), mpdata_V( 
            Cy(idx(ir, jm + grid->p_half, km)), Cy(idx(il, jm + grid->p_half, km)), /* Vru, Vlu */ 
            Cy(idx(ir, jm - grid->m_half, km)), Cy(idx(il, jm - grid->m_half, km))  /* Vrd, Vld */ 
          ) 
        ) + 
        mpdata_CB( 
          (*psi[n])(idx(ir, jm, km+1)), (*psi[n])(idx(il, jm, km+1)), /* pru, plu */ 
          (*psi[n])(idx(ir, jm, km-1)), (*psi[n])(idx(il, jm, km-1)), /* prd, pld */ 
          Cx(idx(ic, jm, km)), mpdata_W( 
            Cz(idx(ir, jm, km + grid->p_half)), Cz(idx(il, jm, km + grid->p_half)), /* Wru, Wlu */ 
            Cz(idx(ir, jm, km - grid->m_half)), Cz(idx(il, jm, km - grid->m_half))  /* Wrd, Wld */ 
          ) 
        )  
      );
      if (third_order) 
      {
        (*C_adf)(adfidx) +=(
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
          ) +
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
          ) -
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
    if (third_order)
    { 
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
    arr<real_t> *psi[], arr<real_t> *[], arr<real_t> *tmp_v[],
    const rng &i, const rng &j, const rng &k, 
    const int n, const int step,
    const arr<real_t> * const Cx, const arr<real_t> * const Cy, const arr<real_t> * const Cz
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

    const arr<real_t> 
      * const Cx_unco = (step < 3 ? Cx : tmp_v[x_old]), 
      *       Cx_corr = (step < 2 ? Cx : tmp_v[x_new]),
      * const Cy_unco = (step < 3 ? Cy : tmp_v[y_old]), 
      *       Cy_corr = (step < 2 ? Cy : tmp_v[y_new]),
      * const Cz_unco = (step < 3 ? Cz : tmp_v[z_old]), 
      *       Cz_corr = (step < 2 ? Cz : tmp_v[z_new]);
     
    *psi[n+1] = *psi[0]; // TODO: at least this should be placed in adv... and the leapfrog & upstream in adv_dimsplit?

    if (true)
    {
      if (step > 1) mpdata_U<idx_ijk>(Cx_corr, psi, n, step, i, j, k, *Cx_unco, *Cy_unco, *Cz_unco);
      adv_upstream<real_t>::template op<idx_ijk>(psi, NULL, NULL, i, j, k, n, 1, Cx_corr, NULL, NULL);
    }
    if (j.first() != j.last())
    {
      if (step > 1) mpdata_U<idx_jki>(Cy_corr, psi, n, step, j, k, i, *Cy_unco, *Cz_unco, *Cx_unco);
      adv_upstream<real_t>::template op<idx_jki>(psi, NULL, NULL, j, k, i, n, 1, Cy_corr, NULL, NULL);
    }
    if (k.first() != k.last())
    {
      if (step > 1) mpdata_U<idx_kij>(Cz_corr, psi, n, step, k, i, j, *Cz_unco, *Cx_unco, *Cy_unco);
      adv_upstream<real_t>::template op<idx_kij>(psi, NULL, NULL, k, i, j, n, 1, Cz_corr, NULL, NULL);
    }
  }
};
#endif
