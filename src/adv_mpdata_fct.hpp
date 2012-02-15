/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011 - January 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @section DESCRIPTION
 *    Flux Corrected Transport (aka non-oscillatory, monotonic, sign-preserving)
 *    option for MPDATA (for the transport of positive scalars only)
 */
#ifndef ADV_MPDATA_FCT_HPP
#  define ADV_MPDATA_FCT_HPP

#  include "adv_mpdata.hpp"

template <typename real_t> 
class adv_mpdata_fct : public adv_mpdata<real_t> 
{
  public: const int stencil_extent() 
  {
    int halo = (adv_mpdata<real_t>::stencil_extent() - 1) / 2; 
    halo += 1; // cf. ii, jj & kk below
    return 1 + 2 * halo;
  }
  public: const int num_vctr_caches() 
  { 
    return 1 + (iord < 3 ? 3 : 6); 
  }
  public: const int num_sclr_caches() { return 2; } 
#  define psi_min (*tmp_s[0])
#  define psi_max (*tmp_s[1])

  private: int iord;
  private: grd_arakawa_c_lorenz<real_t> *grid;
  public: adv_mpdata_fct(grd_arakawa_c_lorenz<real_t> *grid, int iord, bool cross_terms, bool third_order) 
    : adv_mpdata<real_t>(grid, iord, cross_terms, third_order), iord(iord), grid(grid)
  {}

  public: void op3D(
    mtx::arr<real_t> *psi[], mtx::arr<real_t> *tmp_s[], mtx::arr<real_t> *tmp_v[],
    const mtx::idx &ijk,
    const int n, const int step,
    const mtx::arr<real_t> * const Cx, const mtx::arr<real_t> * const Cy, const mtx::arr<real_t> * const Cz
  )
  {
    assert(min((*psi[n])(ijk)) >= real_t(0)); // accepting positive scalars only

    mtx::rng 
      ii = mtx::rng(ijk.lbound(mtx::i) - 1, ijk.ubound(mtx::i) + 1),
      jj = mtx::rng(ijk.lbound(mtx::j) - 1, ijk.ubound(mtx::j) + 1),
      kk = mtx::rng(ijk.lbound(mtx::k) - 1, ijk.ubound(mtx::k) + 1);

    assert(isfinite(sum((*psi[n])(ii  , jj  , kk  ))));
    assert(isfinite(sum((*psi[n])(ii-1, jj  , kk  ))));
    assert(isfinite(sum((*psi[n])(ii+1, jj  , kk  ))));
    assert(isfinite(sum((*psi[n])(ii  , jj-1, kk  ))));
    assert(isfinite(sum((*psi[n])(ii  , jj+1, kk  ))));
    assert(isfinite(sum((*psi[n])(ii  , jj  , kk-1))));
    assert(isfinite(sum((*psi[n])(ii  , jj  , kk+1))));

#  define mpdata_fct_minmax(fun, psi_, n_, i_, j_, k_) fun( \
     (*psi_[n_])(i_  ,j_  ,k_  ), fun( \
     (*psi_[n_])(i_-1,j_  ,k_  ), fun( \
     (*psi_[n_])(i_+1,j_  ,k_  ), fun( \
     (*psi_[n_])(i_  ,j_-1,k_  ), fun( \
     (*psi_[n_])(i_  ,j_+1,k_  ), fun( \
     (*psi_[n_])(i_  ,j_  ,k_-1), \
     (*psi_[n_])(i_  ,j_  ,k_+1)))))) \
   ) 

    ///
    /// \f$ \psi^{max}_{i}=max_{I}(\psi^{n}_{i-1},\psi^{n}_{i},\psi^{n}_{i+1},\psi^{*}_{i-1},\psi^{*}_{i},\psi^{*}_{i+1}) \f$ \n
    /// \f$ \psi^{min}_{i}=min_{I}(\psi^{n}_{i-1},\psi^{n}_{i},\psi^{n}_{i+1},\psi^{*}_{i-1},\psi^{*}_{i},\psi^{*}_{i+1}) \f$ \n
    /// eq.(20a, 20b) in Smolarkiewicz & Grabowski 1990 (J.Comp.Phys.,86,355-375)
    if (step == 1)
    {
      // calculating psi_min and psi_max from the previous time step
      psi_min(ii,jj,kk) = mpdata_fct_minmax(min, psi, n, ii, jj, kk);
      psi_max(ii,jj,kk) = mpdata_fct_minmax(max, psi, n, ii, jj, kk);
      // performing standard upstream advection
      adv_upstream<real_t>::op3D(psi, NULL, NULL, ijk, n, 1, Cx, Cy, Cz);
    }
    else
    {
      // calculating psi_min and psi_max from the previous time step and previous iord
      psi_min(ii,jj,kk) = min(psi_min(ii,jj,kk), mpdata_fct_minmax(min, psi, n, ii, jj, kk));
      psi_max(ii,jj,kk) = max(psi_max(ii,jj,kk), mpdata_fct_minmax(max, psi, n, ii, jj, kk));

      // calculating Cx_mon, Cy_mon, Cz_mon
      int 
        x_old = -1, y_old = -1, z_old = -1, 
        x_new =  1, y_new =  2, z_new =  3; 
      if (iord > 2)
      {   
        if (step > 2)
        {
          x_old = 1, y_old = 2, z_old = 3;
          if (step % 2 == 0)  
            x_old += 3, y_old += 3, z_old += 3;
          else 
            x_new += 3, y_new += 3, z_new += 3;
        }
      }

      const mtx::arr<real_t> 
        * const Cx_unco = (step < 3 ? Cx : tmp_v[x_old]), 
        *       Cx_corr = tmp_v[x_new],
        *       Cx_mono = tmp_v[0],
        * const Cy_unco = (step < 3 ? Cy : tmp_v[y_old]), 
        *       Cy_corr = tmp_v[y_new],
        *       Cy_mono = tmp_v[0],
        * const Cz_unco = (step < 3 ? Cz : tmp_v[z_old]), 
        *       Cz_corr = tmp_v[z_new],
        *       Cz_mono = tmp_v[0];

      this->template mpdata_U<mtx::idx_ijk>(Cx_corr, psi, n, step, ii, jj, kk, *Cx_unco, *Cy_unco, *Cz_unco);
      this->template mpdata_U<mtx::idx_jki>(Cy_corr, psi, n, step, jj, kk, ii, *Cy_unco, *Cz_unco, *Cx_unco);
      this->template mpdata_U<mtx::idx_kij>(Cz_corr, psi, n, step, kk, ii, jj, *Cz_unco, *Cx_unco, *Cy_unco); 
   
      // performing upstream advection using the ''monotonic'' velocities (logic from adv::op3D)
      *psi[n+1] = *psi[0]; // TODO: at least this should be placed in adv... and the leapfrog & upstream in adv_dimsplit?
      if (ijk.i_spans)
      {
        fct_helper<mtx::idx_ijk>(psi, tmp_s, tmp_v, *Cx_mono, *Cx_corr, *Cy_corr, *Cz_corr, ijk.i, ijk.j, ijk.k, n);
        adv_upstream<real_t>::template op<mtx::idx_ijk>(psi, NULL, NULL, ijk.i, ijk.j, ijk.k, n, 1, Cx_mono, NULL, NULL); 
      }
      if (ijk.j_spans)
      {
        fct_helper<mtx::idx_jki>(psi, tmp_s, tmp_v, *Cy_mono, *Cy_corr, *Cz_corr, *Cx_corr, ijk.j, ijk.k, ijk.i, n); 
        adv_upstream<real_t>::template op<mtx::idx_jki>(psi, NULL, NULL, ijk.j, ijk.k, ijk.i, n, 1, Cy_mono, NULL, NULL); 
      }
      if (ijk.k_spans)
      {
        fct_helper<mtx::idx_kij>(psi, tmp_s, tmp_v, *Cz_mono, *Cz_corr, *Cx_corr, *Cy_corr, ijk.k, ijk.i, ijk.j, n); 
        adv_upstream<real_t>::template op<mtx::idx_kij>(psi, NULL, NULL, ijk.k, ijk.i, ijk.j, n, 1, Cz_mono, NULL, NULL); 
      }
    }
  }

  private:
  template <class idx>
  void fct_helper(
    mtx::arr<real_t> *psi[], mtx::arr<real_t> *tmp_s[], mtx::arr<real_t> *[], 
    const mtx::arr<real_t> &C_mon_x, 
    const mtx::arr<real_t> &C_adf_x, const mtx::arr<real_t> &C_adf_y, const mtx::arr<real_t> &C_adf_z, 
    const mtx::rng &i, const mtx::rng &j, const mtx::rng &k, int n
  )  
  {
///
/// \f$ \beta^{\uparrow}_{i} = \frac { \psi^{max}_{i}- \psi^{*}_{i} }
/// { \sum\limits_{I} \frac{\Delta t}{\Delta x^{I}} \left( [u^{I}_{i-1/2}]^{+} \psi^{*}_{i-1} - 
/// [u^{I}_{i+1/2}]^{-} \psi^{*}_{i+1} \right)  } \f$ \n
/// \f$ \beta^{\downarrow}_{i} = \frac { \psi^{*}_{i}- \psi^{min}_{i} }
/// { \sum\limits_{I} \frac{\Delta t}{\Delta x^{I}} \left( [u^{I}_{i+1/2}]^{+} \psi^{*}_{i} - 
/// [u^{I}_{i-1/2}]^{-} \psi^{*}_{i} \right)  } \f$ \n
/// eq.(19a, 19b) in Smolarkiewicz & Grabowski 1990 (J.Comp.Phys.,86,355-375)

#  define mpdata_fct_beta_up(_i, _j, _k) this->mpdata_frac( \
     psi_max(idx(_i,_j,_k)) - (*psi[n])(idx(_i,_j,_k)), \
     mpdata_pospart(C_adf_x(idx(_i - grid->m_half,_j,_k))) * (*psi[n])(idx(_i-1,_j,_k)) - \
     mpdata_negpart(C_adf_x(idx(_i + grid->p_half,_j,_k))) * (*psi[n])(idx(_i+1,_j,_k)) + \
     mpdata_pospart(C_adf_y(idx(_i,_j - grid->m_half,_k))) * (*psi[n])(idx(_i,_j-1,_k)) - \
     mpdata_negpart(C_adf_y(idx(_i,_j + grid->p_half,_k))) * (*psi[n])(idx(_i,_j+1,_k)) + \
     mpdata_pospart(C_adf_z(idx(_i,_j,_k - grid->m_half))) * (*psi[n])(idx(_i,_j,_k-1)) - \
     mpdata_negpart(C_adf_z(idx(_i,_j,_k + grid->p_half))) * (*psi[n])(idx(_i,_j,_k+1))   \
   )
#  define mpdata_fct_beta_dn(_i, _j, _k) this->mpdata_frac(\
     (*psi[n])(idx(_i,_j,_k)) - psi_min(idx(_i,_j,_k)), \
     mpdata_pospart(C_adf_x(idx(_i + grid->p_half,_j,_k))) * (*psi[n])(idx(_i,_j,_k)) - \
     mpdata_negpart(C_adf_x(idx(_i - grid->m_half,_j,_k))) * (*psi[n])(idx(_i,_j,_k)) + \
     mpdata_pospart(C_adf_y(idx(_i,_j + grid->p_half,_k))) * (*psi[n])(idx(_i,_j,_k)) - \
     mpdata_negpart(C_adf_y(idx(_i,_j - grid->m_half,_k))) * (*psi[n])(idx(_i,_j,_k)) + \
     mpdata_pospart(C_adf_z(idx(_i,_j,_k + grid->p_half))) * (*psi[n])(idx(_i,_j,_k)) - \
     mpdata_negpart(C_adf_z(idx(_i,_j,_k - grid->m_half))) * (*psi[n])(idx(_i,_j,_k))   \
   )
    /// nonoscillatory antidiffusive velocity: \n
    /// \f$ U^{MON}_{i+1/2}=min(1,\beta ^{\downarrow}_i,\beta ^{\uparrow} _{i+1})[U_{i+1/2}]^{+} 
    /// + min(1,\beta^{\uparrow}_{i},\beta^{\downarrow}_{i+1/2})[u_{i+1/2}]^{-} \f$ \n
    /// where \f$ [\cdot]^{+}=max(\cdot,0) \f$ and \f$ [\cdot]^{-}=min(\cdot,0) \f$ \n
    /// eq.(18) in Smolarkiewicz & Grabowski 1990 (J.Comp.Phys.,86,355-375)

    // as in mpdata_U, we compute u_{i+1/2} for iv=(i-1, ... i) instead of u_{i+1/2} and u_{i-1/2} for all i
    mtx::rng iv(i.first()-1, i.last());

    assert(isfinite(sum((*psi[n])(idx(iv    ,j ,  k  )))));
    assert(isfinite(sum((*psi[n])(idx(iv+1  ,j ,  k  )))));
    assert(isfinite(sum((*psi[n])(idx(iv  -1,j ,  k  )))));
    assert(isfinite(sum((*psi[n])(idx(iv+1+1,j ,  k  )))));
    assert(isfinite(sum((*psi[n])(idx(iv,    j-1, k  )))));
    assert(isfinite(sum((*psi[n])(idx(iv+1,  j+1, k  )))));
    assert(isfinite(sum((*psi[n])(idx(iv,    j  , k-1)))));
    assert(isfinite(sum((*psi[n])(idx(iv+1,  j  , k+1)))));

    assert(isfinite(sum(C_adf_x(idx(iv     + grid->p_half,j,k)))));
    assert(isfinite(sum(C_adf_x(idx(iv     - grid->m_half,j,k)))));
    assert(isfinite(sum(C_adf_x(idx(iv + 1 + grid->p_half,j,k)))));
    assert(isfinite(sum(C_adf_x(idx(iv + 1 - grid->m_half,j,k)))));

    assert(isfinite(sum(C_adf_y(idx(iv,     j + grid->p_half,k)))));
    assert(isfinite(sum(C_adf_y(idx(iv,     j - grid->m_half,k)))));
    assert(isfinite(sum(C_adf_y(idx(iv + 1, j + grid->p_half,k)))));
    assert(isfinite(sum(C_adf_y(idx(iv + 1, j - grid->m_half,k)))));

    assert(isfinite(sum(C_adf_z(idx(iv,     j, k + grid->p_half)))));
    assert(isfinite(sum(C_adf_z(idx(iv,     j, k - grid->m_half)))));
    assert(isfinite(sum(C_adf_z(idx(iv + 1, j, k + grid->p_half)))));
    assert(isfinite(sum(C_adf_z(idx(iv + 1, j, k - grid->m_half)))));

    C_mon_x(idx(iv + grid->p_half,j,k)) = C_adf_x(idx(iv + grid->p_half,j,k)) * where(
      C_adf_x(idx(iv + grid->p_half,j,k)) > 0,
      min(1, min(mpdata_fct_beta_dn(iv, j, k), mpdata_fct_beta_up(iv+1, j, k))),
      min(1, min(mpdata_fct_beta_up(iv, j, k), mpdata_fct_beta_dn(iv+1, j, k)))
    );
  }
};
#endif
