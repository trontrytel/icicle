/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011 - March 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @brief definition of the adv_mpdata class with an implementation of the MPDATA advection scheme
 */
#ifndef ADV_MPDATA_HPP
#  define ADV_MPDATA_HPP

#  include "adv_upstream.hpp"
#  include "grd_arakawa-c-lorenz.hpp"

/** @brief
 *  C++ implementation of the second- and third-order accurate
 *  MPDATA scheme for solenoidal flows of scalar fields on a carthesian
 *  1D, 2D and 3D Arakawa-C grid 
 */
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
    if (iord < 3 && third_order) warning_macro("third-order accuracy needs iord >= 3")
  }

  public: 
  template <class aon_class>
  class op3D : public aon_class, public adv<real_t>::op3D
  {
    // nested class for storing indices
    protected:
    template <class idx>
    class indices 
    {   
      public: bool do_x, do_y, do_z;
      private: int iord_halo_yz, iord_halo_x;
      private: mtx::rng im, jm, km; // instead of computing u_{i+1/2} and u_{i-1/2} for all i we compute u_{i+1/2} for im=(i-1, ... i)
      private: mtx::rng ir, ic, il; // forward-in-space perspective
      public: idx adfidx, ir_jm_km, il_jmm1_km, ic_jm_km, il_jm_km, imm1_jm_km, im_jm_km, imp1_jm_km, imp2_jm_km,
        il_jm_kmmh, ir_jm_kmmh, il_jm_kmph, ir_jm_kmph, il_jmmh_km, ir_jmmh_km, il_jmph_km, ir_jmph_km, imp1_jmm1_kmm1,
        imp1_jmp1_kmm1, im_jmm1_kmm1, im_jmp1_kmm1, imp1_jmm1_kmp1, imp1_jmp1_kmp1, im_jmm1_kmp1, im_jmp1_kmp1,
        imp1_jm_kmm1, imp1_jm_kmp1, im_jm_kmm1, im_jm_kmp1, imp1_jmm1_km, imp1_jmp1_km, im_jmm1_km, im_jmp1_km,
        il_jm_kmm1, ir_jm_kmm1, il_jm_kmp1, ir_jm_kmp1, ir_jmm1_km, il_jmp1_km, ir_jmp1_km;
      public: indices(
        const mtx::rng &i, 
        const mtx::rng &j, 
        const mtx::rng &k, 
        const grd_arakawa_c_lorenz<real_t> &grid,
        int cross_terms,
        int iord, int step
      ) :
        do_x(i.first() != i.last()), 
        do_y(j.first() != j.last()),
        do_z(k.first() != k.last()),
        iord_halo_yz((iord > 2 && cross_terms) ? iord - step : 0),
        iord_halo_x((iord > 3 && cross_terms && iord != step) ? iord - step - 1 : 0),
        im(i.first() - 1 - iord_halo_x, i.last() + iord_halo_x),
        jm(j.first() - iord_halo_yz, j.last() + iord_halo_yz), 
        km(k.first() - iord_halo_yz, k.last() + iord_halo_yz),
        ir(im + 1),           // right
        ic(im + grid.p_half), // center
        il(im),               // left
        adfidx(
          mtx::rng(i.first() - grid.m_half - iord_halo_x, i.last() + grid.p_half + iord_halo_x),
          jm, 
          km
        ),
        ir_jm_km(ir, jm, km),
        il_jmm1_km(il, jm-1, km),
        ic_jm_km(ic, jm, km),
        il_jm_km(il, jm, km),
        imm1_jm_km(im-1, jm, km),
        im_jm_km(im, jm, km),
        imp1_jm_km(im+1, jm, km),
        imp2_jm_km(im+2, jm, km),
        il_jm_kmmh(il, jm, km - grid.m_half),
        ir_jm_kmmh(ir, jm, km - grid.m_half),
        il_jm_kmph(il, jm, km + grid.p_half),
        ir_jm_kmph(ir, jm, km + grid.p_half),
        il_jmmh_km(il, jm - grid.m_half, km),
        ir_jmmh_km(ir, jm - grid.m_half, km),
        il_jmph_km(il, jm + grid.p_half, km),
        ir_jmph_km(ir, jm + grid.p_half, km),
        imp1_jmm1_kmm1(im + 1, jm - 1, km - 1),
        imp1_jmp1_kmm1(im + 1, jm + 1, km - 1),
        im_jmm1_kmm1(im, jm - 1, km - 1),
        im_jmp1_kmm1(im, jm + 1, km - 1),
        imp1_jmm1_kmp1(im + 1, jm - 1, km + 1),
        imp1_jmp1_kmp1(im + 1, jm + 1, km + 1),
        im_jmm1_kmp1(im, jm - 1, km + 1),
        im_jmp1_kmp1(im, jm + 1, km + 1), 
        imp1_jm_kmm1(im + 1, jm, km - 1),
        imp1_jm_kmp1(im + 1, jm, km + 1),
        im_jm_kmm1(im, jm, km - 1), 
        im_jm_kmp1(im, jm, km + 1), 
        imp1_jmm1_km(im + 1, jm - 1, km), 
        imp1_jmp1_km(im + 1, jm + 1, km), 
        im_jmm1_km(im, jm - 1, km), 
        im_jmp1_km(im, jm + 1, km),
        il_jm_kmm1(il, jm, km - 1), 
        ir_jm_kmm1(ir, jm, km - 1), 
        il_jm_kmp1(il, jm, km + 1), 
        ir_jm_kmp1(ir, jm, km + 1), 
        ir_jmm1_km(ir, jm - 1, km), 
        il_jmp1_km(il, jm + 1, km),
        ir_jmp1_km(ir, jm + 1, km)
      { } 
    };  

    // members fields
    private: typename adv_upstream<real_t>::op3D upstream;
    private: mtx::arr<real_t> **tmp_s, **tmp_v;
    private: int iord, cross_terms, third_order;
    public: ptr_vector<indices<mtx::idx_ijk>> indcs_x;
    public: ptr_vector<indices<mtx::idx_jki>> indcs_y;
    public: ptr_vector<indices<mtx::idx_kij>> indcs_z;

    // ctor
    public: op3D(
      const mtx::idx &ijk, 
      const grd_arakawa_c_lorenz<real_t> &grid, 
      mtx::arr<real_t> **tmp_s,
      mtx::arr<real_t> **tmp_v,
      int cross_terms, int iord, int third_order
    ) :
      adv<real_t>::op3D(ijk),
      upstream(ijk, grid),
      tmp_s(tmp_s), tmp_v(tmp_v), iord(iord), cross_terms(cross_terms), third_order(third_order)
    { 
      for (int step = 1; step <= iord; ++step)
      {
        indcs_x.push_back(new indices<mtx::idx_ijk>(ijk.i, ijk.j, ijk.k, grid, cross_terms, iord, step));
        indcs_y.push_back(new indices<mtx::idx_jki>(ijk.j, ijk.k, ijk.i, grid, cross_terms, iord, step)); 
        indcs_z.push_back(new indices<mtx::idx_kij>(ijk.k, ijk.i, ijk.j, grid, cross_terms, iord, step));
      }
    } 

    /// multidimensional antidiffusive velocity: 
    /// \f$ 
    ///   \tilde{U}^{I}_{i+1/2} = \left[ 
    ///     |U^{I}_{i+1/2}| \Delta x^{I} - \Delta t (u^{I}_{i+1/2})^2 
    ///   \right] \cdot
    ///   \frac{
    ///     \psi^{*}_{i+1}-\psi^{*}_{i}
    ///   }{
    ///     (\psi^{*}_{i+1}+\psi^{*}_{i}) \Delta x^{I}
    ///   }  
    ///   - \sum\limits_{J=1,J \ne I} 0.5 \Delta t U^{I}_{i+1/2} \bar{U}^{J}_{i+1/2} \cdot
    ///   \frac{
    ///     \psi^{*}_{i+1,j+1}+\psi^{*}_{i,j+1}-\psi^{*}_{i+1,j-1}-\psi^{*}_{i,j-1}
    ///   }{
    ///     \psi^{*}_{i+1,j+1}+\psi^{*}_{i,j+1}+\psi^{*}_{i+1,j-1}+\psi^{*}_{i,j-1}
    ///   } 
    /// \f$
    /// eq. (13-14) in Smolarkiewicz 1984 (J. Comp. Phys.,54,352-362) \n

    public:
    template <class indices>
    void mpdata_U(
      const mtx::arr<real_t> * C_adf,
      const mtx::arr<real_t> * const psi[], const int n, const int step,
      const ptr_vector<indices> &indcs,
      const mtx::arr<real_t> &Cx, 
      const mtx::arr<real_t> &Cy, 
      const mtx::arr<real_t> &Cz
    )
    {
      const indices &idx = indcs[step-1];

      (*C_adf)(idx.adfidx) = (
        adv_mpdata<real_t>::f_CA( 
          this->aon((*psi[n])(idx.ir_jm_km)), 
          this->aon((*psi[n])(idx.il_jm_km)), 
          Cx(idx.ic_jm_km) 
        )
      );  
      if (cross_terms) 
      {
        if (idx.do_y)
        {
          (*C_adf)(idx.adfidx) -= (
            adv_mpdata<real_t>::f_CB( 
              this->aon((*psi[n])(idx.ir_jmp1_km)), 
              this->aon((*psi[n])(idx.il_jmp1_km)), 
              this->aon((*psi[n])(idx.ir_jmm1_km)), 
              this->aon((*psi[n])(idx.il_jmm1_km)), 
              Cx(idx.ic_jm_km), adv_mpdata<real_t>::f_V( 
                Cy(idx.ir_jmph_km), Cy(idx.il_jmph_km), 
                Cy(idx.ir_jmmh_km), Cy(idx.il_jmmh_km)  
              ) 
            )
          );
        }
        if (idx.do_z)
        {
          (*C_adf)(idx.adfidx) -= ( // otherwise Cz is uninitialised!
            adv_mpdata<real_t>::f_CB( 
              this->aon((*psi[n])(idx.ir_jm_kmp1)), 
              this->aon((*psi[n])(idx.il_jm_kmp1)),
              this->aon((*psi[n])(idx.ir_jm_kmm1)), 
              this->aon((*psi[n])(idx.il_jm_kmm1)),
              Cx(idx.ic_jm_km), adv_mpdata<real_t>::f_W( 
                Cz(idx.ir_jm_kmph), Cz(idx.il_jm_kmph), 
                Cz(idx.ir_jm_kmmh), Cz(idx.il_jm_kmmh)  
              ) 
            )  
          );
        }
        if (third_order) 
        {
          if (idx.do_y)
          {
            (*C_adf)(idx.adfidx) += ( // otherwise Cy is uninitialised
              adv_mpdata<real_t>::f_3rd_xy(
                this->aon((*psi[n])(idx.im_jmp1_km)), 
                this->aon((*psi[n])(idx.im_jmm1_km)), 
                this->aon((*psi[n])(idx.imp1_jmp1_km)), 
                this->aon((*psi[n])(idx.imp1_jmm1_km)), 
                Cx(idx.ic_jm_km), // U, 
                adv_mpdata<real_t>::f_V( // V
                  Cy(idx.ir_jmph_km), Cy(idx.il_jmph_km), 
                  Cy(idx.ir_jmmh_km), Cy(idx.il_jmmh_km)  
                ) 
              )
            );
          }
          if (idx.do_z)
          {
            (*C_adf)(idx.adfidx) += ( // otherwise Cz is uninitialised
              adv_mpdata<real_t>::f_3rd_xz(
                this->aon((*psi[n])(idx.im_jm_kmp1)), 
                this->aon((*psi[n])(idx.im_jm_kmm1)), 
                this->aon((*psi[n])(idx.imp1_jm_kmp1)), 
                this->aon((*psi[n])(idx.imp1_jm_kmm1)), 
                Cx(idx.ic_jm_km), // U, 
                adv_mpdata<real_t>::f_W( // W
                  Cz(idx.ir_jm_kmph), Cz(idx.il_jm_kmph), /* Wru, Wlu */ 
                  Cz(idx.ir_jm_kmmh), Cz(idx.il_jm_kmmh)  /* Wrd, Wld */ 
                ) 
              )
            ); 
          }
          if (idx.do_y && idx.do_z)
          {
            (*C_adf)(idx.adfidx) -= ( // otherwise Cx & Cz are uninitialised
              adv_mpdata<real_t>::f_3rd_yz(
                this->aon((*psi[n])(idx.im_jmp1_kmp1)), 
                this->aon((*psi[n])(idx.im_jmm1_kmp1)),
                this->aon((*psi[n])(idx.imp1_jmp1_kmp1)), 
                this->aon((*psi[n])(idx.imp1_jmm1_kmp1)), 
                this->aon((*psi[n])(idx.im_jmp1_kmm1)), 
                this->aon((*psi[n])(idx.im_jmm1_kmm1)), 
                this->aon((*psi[n])(idx.imp1_jmp1_kmm1)), 
                this->aon((*psi[n])(idx.imp1_jmm1_kmm1)), 
                Cx(idx.ic_jm_km), // U
                adv_mpdata<real_t>::f_V( // V
                  Cy(idx.ir_jmph_km), Cy(idx.il_jmph_km), /* Vru, Vlu */ 
                  Cy(idx.ir_jmmh_km), Cy(idx.il_jmmh_km)  /* Vrd, Vld */ 
                ), 
                adv_mpdata<real_t>::f_W( // W
                  Cz(idx.ir_jm_kmph), Cz(idx.il_jm_kmph), /* Wru, Wlu */ 
                  Cz(idx.ir_jm_kmmh), Cz(idx.il_jm_kmmh)  /* Wrd, Wld */ 
                ) 
              )
            );
          }
        }
      }
      if (third_order)
      { 
        (*C_adf)(idx.adfidx) +=(
          adv_mpdata<real_t>::f_3rd_xx(
            this->aon((*psi[n])(idx.imp2_jm_km)), 
            this->aon((*psi[n])(idx.imp1_jm_km)), 
            this->aon((*psi[n])(idx.im_jm_km)), 
            this->aon((*psi[n])(idx.imm1_jm_km)), 
            Cx(idx.ic_jm_km)
          )
        );  
      }
    }

    /// \f$ 
    ///   \psi^{n+1}_{i} = \psi^{*}_{i} -\sum\limits_{I} \left[ 
    ///     F^{I}( \psi^{*}_{i}, \psi^{*}_{i+1}, \tilde{U}^{I}_{i+1/2} ) -
    ///     F^{I}( \psi^{*}_{i-1}, \psi^{*}_{i}, \tilde{U}^{I}_{i-1/2} ) 
    ///   \right] 
    /// \f$
    /// where \f$ I \f$ denotes the sum over all dimensions and \f$ \tilde{U} \f$ is 
    /// the multidimensional antidiffusive velocity eq. (12) in Smolarkiewicz 1984 (J. Comp. Phys.,54,352-362) 
    public: void operator()(
      mtx::arr<real_t> *psi[],
      const int n,
      const int step,
      const mtx::arr<real_t> * const Cx, 
      const mtx::arr<real_t> * const Cy, 
      const mtx::arr<real_t> * const Cz
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
      if (this->do_x())
      {
        if (step > 1) mpdata_U(Cx_corr, psi, n, step, indcs_x, *Cx_unco, *Cy_unco, *Cz_unco);
        upstream.op1D(psi, upstream.indcs_x, n, 1, Cx_corr, NULL, NULL);
      }
      if (this->do_y())
      {
        if (step > 1) mpdata_U(Cy_corr, psi, n, step, indcs_y, *Cy_unco, *Cz_unco, *Cx_unco);
        upstream.op1D(psi, upstream.indcs_y, n, 1, Cy_corr, NULL, NULL);
      }
      if (this->do_z())
      {
        if (step > 1) mpdata_U(Cz_corr, psi, n, step, indcs_z, *Cz_unco, *Cx_unco, *Cy_unco);
        upstream.op1D(psi, upstream.indcs_z, n, 1, Cz_corr, NULL, NULL);
      }
    }
  };

  // on-demand support for transporting fields of variable sign
  // see section 3.2, subsection (4) in PKS & LGM 1998 
  protected: class aon_nil 
  { 
    protected: template<class T> T aon(const T &x) 
    { 
      return x; 
    } 
  };
  protected: class aon_abs 
  { 
    protected: template<class T> auto aon(const T &x) -> decltype(abs(x))
    { 
      return abs(x); 
    } 
  };

  public: typename adv<real_t>::op3D *factory(
    const mtx::idx &ijk,
    mtx::arr<real_t> **tmp_s,
    mtx::arr<real_t> **tmp_v,
    bool positive_definite
  )
  {
    if (positive_definite)
      return new op3D<aon_nil>(ijk, *grid, tmp_s, tmp_v, cross_terms, iord, third_order);
    else
      return new op3D<aon_abs>(ijk, *grid, tmp_s, tmp_v, cross_terms, iord, third_order);
  }
    // TODO: make it an option for the constructor (and recode with functors)
#    ifdef MPDATA_FRAC_EPSILON
    protected: mtx_expr_2arg_macro(frac, num, den, 
      num / (den + mtx::eps<real_t>())
    )
#    else
    protected: mtx_expr_2arg_macro(frac, num, den, 
      where(den > real_t(0), num / den, real_t(0))
    ) 
#    endif

    // macros for 2nd order terms:
    private: mtx_expr_2arg_macro(f_A, pr, pl,
      adv_mpdata<real_t>::frac(pr - pl, pr + pl)
    ) 

    private: mtx_expr_4arg_macro(f_B, pru, plu, prd, pld,
      real_t(.5) * adv_mpdata<real_t>::frac(pru + plu - prd - pld, pru + plu + prd + pld)
    )

    private: mtx_expr_4arg_macro(f_V, Vru, Vlu, Vrd, Vld, 
      real_t(.25) * (Vru + Vlu + Vrd + Vld)
    )

    private: mtx_expr_4arg_macro(f_W, Wru, Wlu, Wrd, Wld,
      adv_mpdata<real_t>::f_V(Wru, Wlu, Wrd, Wld)
    )

    private: mtx_expr_3arg_macro(f_CA, pr, pl, U,
      (abs(U) - pow(U,2)) * adv_mpdata<real_t>::f_A(pr, pl)
    )

    private: mtx_expr_6arg_macro(f_CB, pru, plu, prd, pld, U, V,
      U * V * adv_mpdata<real_t>::f_B(pru, plu, prd, pld)
    )

    // macros for 3rd order terms:
    /// first term from eq. (36) from Smolarkiewicz & Margolin 1998 (with G=1)
    /// \f$ 
    ///   \frac{(\delta x)^2}{6} \left( 3U|U| - 2U^3 - U \right) 
    ///   \frac{1}{\psi} \frac{\partial^2 \psi}{\partial x^2} \approx
    ///   (\ldots) \cdot \frac{1}{3} \cdot 
    ///   \frac{\psi_{i+2} - \psi_{i+1} - \psi{i} + \psi{i-1}}{\psi_{i+2} + \psi_{i+1} + \psi{i} + \psi{i-1}}
    /// \f$ \n 
    private: mtx_expr_5arg_macro(f_3rd_xx, pp2, pp1, pp0, pm1, U, 
      (real_t(3) * U * abs(U) - real_t(2) * pow(U,3) - U) 
      / real_t(3) 
      * adv_mpdata<real_t>::frac(pp2 - pp1 - pp0 + pm1, pp2 + pp1 + pp0 + pm1) 
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
    private: mtx_expr_6arg_macro(f_3rd_xy, pip0jp1, pip0jm1, pip1jp1, pip1jm1, U, V,
      V * (abs(U) - real_t(2) * pow(U,2)) 
      * adv_mpdata<real_t>::frac( 
        pip1jp1 - pip0jp1 - pip1jm1 + pip0jm1, 
        pip1jp1 + pip0jp1 + pip1jm1 + pip0jm1 
      ) 
    )

    /// third term from eq. (36) from Smolarkiewicz & Margolin 1998 (with G=1)
    private: mtx_expr_6arg_macro(f_3rd_xz, pip0kp1, pip0km1, pip1kp1, pip1km1, U, W,
      adv_mpdata<real_t>::f_3rd_xy(pip0kp1, pip0km1, pip1kp1, pip1km1, U, W)
    )

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
    private: mtx_expr_11arg_macro(f_3rd_yz,
      pip0jp1kp1, pip0jm1kp1, pip1jp1kp1, pip1jm1kp1, 
      pip0jp1km1, pip0jm1km1, pip1jp1km1, pip1jm1km1, 
      U, V, W, 
      real_t(2./3.) * U * V * W * adv_mpdata<real_t>::frac( 
        pip0jp1kp1 - pip0jm1kp1 + pip1jp1kp1 - pip1jm1kp1 - 
        pip0jp1km1 + pip0jm1km1 - pip1jp1km1 + pip1jm1km1, 
        pip0jp1kp1 + pip0jm1kp1 + pip1jp1kp1 + pip1jm1kp1 + 
        pip0jp1km1 + pip0jm1km1 + pip1jp1km1 + pip1jm1km1 
      ) 
    )

};
#endif
