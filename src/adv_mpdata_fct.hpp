/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 *  @section DESCRIPTION
 *    Flux Corrected Transport aka non-oscillatory aka monotonic
 *    option for MPDATA (for the transport of positive scalars only)
 */
#ifndef ADV_MPDATA_FCT_HPP
#  define ADV_MPDATA_FCT_HPP

#  include "adv_mpdata.hpp"

template <typename real_t> 
class adv_mpdata_fct : public adv_mpdata<real_t> 
{
  public: const int stencil_extent() { return 5; }
  public: const int num_vctr_caches() { return 6; }
#  define C_adf(dim) (*tmp_v[0+dim])
#  define C_mon(dim) (*tmp_v[3+dim])
  public: const int num_sclr_caches() { return 2; }
#  define psi_min (*tmp_s[0])
#  define psi_max (*tmp_s[1])

  private: grd_arakawa_c_lorenz<real_t> *grid;
  public: adv_mpdata_fct(grd_arakawa_c_lorenz<real_t> *grid, int iord) 
    : adv_mpdata<real_t>(grid, iord), grid(grid)
  {}

  public: void op3D(
    arr<real_t> *psi[], arr<real_t> *tmp_s[], arr<real_t> *tmp_v[],
    const Range &i, const Range &j, const Range &k,
    const int n, const int step,
    arr<real_t> &Cx, arr<real_t> &Cy, arr<real_t> &Cz
  )
  {
    assert(min(*psi[n]) >= real_t(0)); // accepting positive scalars only

    Range ii = Range(i.first() - 1, i.last() + 1),
          jj = Range(j.first() - 1, j.last() + 1),
          kk = Range(k.first() - 1, k.last() + 1);
#  define mpdata_fct_minmax(fun, psi_, n_, i_, j_, k_) ::fun( \
     (*psi_[n_])(i_  ,j_  ,k_  ), ::fun( \
     (*psi_[n_])(i_-1,j_  ,k_  ), ::fun( \
     (*psi_[n_])(i_+1,j_  ,k_  ), ::fun( \
     (*psi_[n_])(i_  ,j_-1,k_  ), ::fun( \
     (*psi_[n_])(i_  ,j_+1,k_  ), ::fun( \
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
      adv_upstream<real_t>::op3D(psi, NULL, NULL, i, j, k, n, 1, Cx, Cy, Cz);
    }
    else
    {
      // calculating psi_min and psi_max from the previous time step and previous iord
      psi_min(ii,jj,kk) = ::min(psi_min(ii,jj,kk), mpdata_fct_minmax(min, psi, n, ii, jj, kk));
      psi_max(ii,jj,kk) = ::max(psi_max(ii,jj,kk), mpdata_fct_minmax(max, psi, n, ii, jj, kk));
      // calculating Cx_mon, Cy_mon, Cz_mon
      this->template mpdata_U<idx_ijk>(&C_adf(0), psi, n, i, j, k, Cx, Cy, Cz);
      this->template mpdata_U<idx_jki>(&C_adf(1), psi, n, j, k, i, Cy, Cz, Cx); // TODO only if 3D?
      this->template mpdata_U<idx_kij>(&C_adf(2), psi, n, k, i, j, Cz, Cx, Cy); // TODO only for 2D and 3D?
      // TODO: this could be implemented with just one cache for C_mon!
      fct_helper<idx_ijk>(psi, tmp_s, tmp_v, C_mon(0), C_adf(0), C_adf(1), C_adf(2), i, j, k, n);
      fct_helper<idx_jki>(psi, tmp_s, tmp_v, C_mon(1), C_adf(1), C_adf(2), C_adf(0), j, k, i, n); // TODO: only if 3D?
      fct_helper<idx_kij>(psi, tmp_s, tmp_v, C_mon(2), C_adf(2), C_adf(0), C_adf(1), k, i, j, n); // TODO: only for 2D and 3D?
      // performing upstream advection using the calculated ''monotonic'' velocities
      adv_upstream<real_t>::op3D(psi, NULL, NULL, i, j, k, n, 1, C_mon(0), C_mon(1), C_mon(2));
    }
  }

  private:
  template <class idx>
  void fct_helper(
    arr<real_t> *psi[], arr<real_t> *tmp_s[], arr<real_t> *tmp_v[], 
    const arr<real_t> &C_mon_x, 
    const arr<real_t> &C_adf_x, const arr<real_t> &C_adf_y, const arr<real_t> &C_adf_z, 
    const Range &i, const Range &j, const Range &k, int n
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

#  define mpdata_fct_beta_up(_i, _j, _k) mpdata_frac( \
     psi_max(idx(_i,_j,_k)) - (*psi[n])(idx(_i,_j,_k)), \
     (::max(real_t(0),C_adf_x(idx(_i - grid->m_half,_j,_k)))) * (*psi[n])(idx(_i-1,_j,_k)) - \
     (::min(real_t(0),C_adf_x(idx(_i + grid->p_half,_j,_k)))) * (*psi[n])(idx(_i+1,_j,_k)) + \
     (::max(real_t(0),C_adf_y(idx(_i,_j - grid->m_half,_k)))) * (*psi[n])(idx(_i,_j-1,_k)) - \
     (::min(real_t(0),C_adf_y(idx(_i,_j + grid->p_half,_k)))) * (*psi[n])(idx(_i,_j+1,_k)) + \
     (::max(real_t(0),C_adf_z(idx(_i,_j,_k - grid->m_half)))) * (*psi[n])(idx(_i,_j,_k-1)) - \
     (::min(real_t(0),C_adf_z(idx(_i,_j,_k + grid->p_half)))) * (*psi[n])(idx(_i,_j,_k+1))   \
   )
#  define mpdata_fct_beta_dn(_i, _j, _k) mpdata_frac(\
     (*psi[n])(idx(_i,_j,_k)) - psi_min(idx(_i,_j,_k)), \
     (::max(real_t(0),C_adf_x(idx(_i + grid->p_half,_j,_k)))) * (*psi[n])(idx(_i,_j,_k)) - \
     (::min(real_t(0),C_adf_x(idx(_i - grid->m_half,_j,_k)))) * (*psi[n])(idx(_i,_j,_k)) + \
     (::max(real_t(0),C_adf_y(idx(_i,_j + grid->p_half,_k)))) * (*psi[n])(idx(_i,_j,_k)) - \
     (::min(real_t(0),C_adf_y(idx(_i,_j - grid->m_half,_k)))) * (*psi[n])(idx(_i,_j,_k)) + \
     (::max(real_t(0),C_adf_z(idx(_i,_j,_k + grid->p_half)))) * (*psi[n])(idx(_i,_j,_k)) - \
     (::min(real_t(0),C_adf_z(idx(_i,_j,_k - grid->m_half)))) * (*psi[n])(idx(_i,_j,_k))   \
   )
    // as in mpdata_U, we compute u_{i+1/2} for iv=(i-1, ... i) instead of u_{i+1/2} and u_{i-1/2} for all i
    Range iv(i.first()-1, i.last());
    /// nonoscillatory antidiffusive velocity: \n
    /// \f$ U^{MON}_{i+1/2}=min(1,\beta ^{\downarrow}_i,\beta ^{\uparrow} _{i+1})[U_{i+1/2}]^{+} 
    /// + min(1,\beta^{\uparrow}_{i},\beta^{\downarrow}_{i+1/2})[u_{i+1/2}]^{-} \f$ \n
    /// where \f$ [\cdot]^{+}=max(\cdot,0) \f$ and \f$ [\cdot]^{-}=min(\cdot,0) \f$ \n
    /// eq.(18) in Smolarkiewicz & Grabowski 1990 (J.Comp.Phys.,86,355-375)

    C_mon_x(idx(iv + grid->p_half,j,k)) = C_adf_x(idx(iv + grid->p_half,j,k)) * where(
      C_adf_x(idx(iv + grid->p_half,j,k)) > 0,
      ::min(1, ::min(mpdata_fct_beta_dn(iv, j, k), mpdata_fct_beta_up(iv+1, j, k))),
      ::min(1, ::min(mpdata_fct_beta_up(iv, j, k), mpdata_fct_beta_dn(iv+1, j, k)))
    );
  }
};
#endif
