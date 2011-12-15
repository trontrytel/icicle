/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 *  @section DESCRIPTION
 *    Flux Corrected Transport aka non-oscillatory aka monotonic
 *    option for MPDATA
 *    only for the transport of positive scalars
 */
#ifndef ADV_MPDATA_FCT_HPP
#  define ADV_MPDATA_FCT_HPP

// TODO: perhaps add an assert() to check if psi has positive values

#  include "adv_mpdata.hpp"

template <typename real_t> 
class adv_mpdata_fct : public adv_mpdata<real_t> 
{
  public: const int stencil_extent() { return 5; }

  public: const int num_vctr_caches() 
  { 
    return adv_mpdata<real_t>::num_vctr_caches() + 3; 
  }
#  define C_mon(dim) (*tmp_v[3+dim])
#  define C_adf(dim) (*tmp_v[0+dim])

  public: const int num_sclr_caches() 
  { 
    return adv_mpdata<real_t>::num_sclr_caches() + 2; 
  }
#  define psi_min (*tmp_s[0])
#  define psi_max (*tmp_s[1])

  private: grd_arakawa_c_lorenz<real_t> *grid;
  public: adv_mpdata_fct(grd_arakawa_c_lorenz<real_t> *grid, int iord) 
    : adv_mpdata<real_t>(grid, iord, true), grid(grid)
  { }

  public: void op3D(
    Array<real_t, 3> *psi[],
    Array<real_t, 3> *tmp_s[],
    Array<real_t, 3> *tmp_v[],
    const Range &i, const Range &j, const Range &k,
    const int n, const int step,
    Array<real_t,3> &Cx, Array<real_t,3> &Cy, Array<real_t,3> &Cz
  )
  {
    adv_mpdata<real_t>::op3D(psi, tmp_s, tmp_v, i, j, k, n, step, Cx, Cy, Cz);

    Range  
      ii = Range(i.first() - 1, i.last() + 1),
      jj = Range(j.first() - 1, j.last() + 1),
      kk = Range(k.first() - 1, k.last() + 1);

    if (step == 1)
    {
      psi_min(ii,jj,kk) = ::min(
        (*psi[n])(ii  ,jj  ,kk  ), ::min(
        (*psi[n])(ii-1,jj  ,kk  ), ::min(
        (*psi[n])(ii+1,jj  ,kk  ), ::min(
        (*psi[n])(ii  ,jj-1,kk  ), ::min(
        (*psi[n])(ii  ,jj+1,kk  ), ::min(
        (*psi[n])(ii  ,jj  ,kk-1), 
        (*psi[n])(ii  ,jj  ,kk+1)))))));
      psi_max(ii,jj,kk) = ::max(
        (*psi[n])(ii  ,jj  ,kk  ), ::max(
        (*psi[n])(ii-1,jj  ,kk  ), ::max(
        (*psi[n])(ii+1,jj  ,kk  ), ::max(
        (*psi[n])(ii  ,jj-1,kk  ), ::max(
        (*psi[n])(ii  ,jj+1,kk  ), ::max(
        (*psi[n])(ii  ,jj  ,kk-1), 
        (*psi[n])(ii  ,jj  ,kk+1)))))));

    }
    else
    {
      // calculating psi_min and psi_max
      psi_min(ii,jj,kk) = ::min(psi_min(ii,jj,kk), ::min(
        (*psi[n])(ii  ,jj  ,kk  ), ::min(
        (*psi[n])(ii-1,jj  ,kk  ), ::min(
        (*psi[n])(ii+1,jj  ,kk  ), ::min(
        (*psi[n])(ii  ,jj-1,kk  ), ::min(
        (*psi[n])(ii  ,jj+1,kk  ), ::min(
        (*psi[n])(ii  ,jj  ,kk-1), 
        (*psi[n])(ii  ,jj  ,kk+1))))))));
      psi_max(i,j,k) = ::max(psi_max(ii,jj,kk), ::max(
        (*psi[n])(ii  ,jj  ,kk  ), ::max(
        (*psi[n])(ii-1,jj  ,kk  ), ::max(
        (*psi[n])(ii+1,jj  ,kk  ), ::max(
        (*psi[n])(ii  ,jj-1,kk  ), ::max(
        (*psi[n])(ii  ,jj+1,kk  ), ::max(
        (*psi[n])(ii  ,jj  ,kk-1), 
        (*psi[n])(ii  ,jj  ,kk+1))))))));

      // calculating Cx_mon, Cy_mon, Cz_mon
      fct_helper<idx_ijk>(psi, tmp_s, tmp_v, C_mon(0), C_adf(0), C_adf(1), C_adf(2), i, j, k, n);
      fct_helper<idx_jki>(psi, tmp_s, tmp_v, C_mon(1), C_adf(1), C_adf(2), C_adf(0), j, k, i, n);
      fct_helper<idx_kij>(psi, tmp_s, tmp_v, C_mon(2), C_adf(2), C_adf(0), C_adf(1), k, i, j, n);

      // TODO adv_upstream<>::op3D...
      adv_mpdata<real_t>::op3D(psi, tmp_s, tmp_v, i, j, k, n, 1, C_mon(0), C_mon(1), C_mon(2));
    }
  }

  private:
  template <class idx>
  void fct_helper(
    Array<real_t, 3> *psi[], 
    Array<real_t, 3> *tmp_s[], 
    Array<real_t, 3> *tmp_v[], 
    const Array<real_t, 3> &C_mon_x, 
    const Array<real_t, 3> &C_adf_x, 
    const Array<real_t, 3> &C_adf_y, 
    const Array<real_t, 3> &C_adf_z, 
    const Range &i, const Range &j, const Range &k,
    int n
  )  
  {
    // instead of computing u_{i+1/2} and u_{i-1/2} for all i
    // we compute u_{i+1/2} for iv=(i-1, ... i)
    Range iv(i.first()-1, i.last());

#    define mpdata_beta_up(_i, _j, _k) mpdata_frac( \
       psi_max(idx(_i,_j,_k)) - (*psi[n])(idx(_i,_j,_k)), \
       (::max(real_t(0),C_adf_x(idx(_i - grid->m_half,_j,_k)))) * (*psi[n])(idx(_i-1,_j,_k)) - \
       (::min(real_t(0),C_adf_x(idx(_i + grid->p_half,_j,_k)))) * (*psi[n])(idx(_i+1,_j,_k)) + \
       (::max(real_t(0),C_adf_y(idx(_i,_j - grid->m_half,_k)))) * (*psi[n])(idx(_i,_j-1,_k)) - \
       (::min(real_t(0),C_adf_y(idx(_i,_j + grid->p_half,_k)))) * (*psi[n])(idx(_i,_j+1,_k)) + \
       (::max(real_t(0),C_adf_z(idx(_i,_j,_k - grid->m_half)))) * (*psi[n])(idx(_i,_j,_k-1)) - \
       (::min(real_t(0),C_adf_z(idx(_i,_j,_k + grid->p_half)))) * (*psi[n])(idx(_i,_j,_k+1))   \
     )
#    define mpdata_beta_dn(_i, _j, _k) mpdata_frac(\
       (*psi[n])(idx(_i,_j,_k)) - psi_min(idx(_i,_j,_k)), \
       (::max(real_t(0),C_adf_x(idx(_i + grid->p_half,_j,_k)))) * (*psi[n])(idx(_i,_j,_k)) - \
       (::min(real_t(0),C_adf_x(idx(_i - grid->m_half,_j,_k)))) * (*psi[n])(idx(_i,_j,_k)) + \
       (::max(real_t(0),C_adf_y(idx(_i,_j + grid->p_half,_k)))) * (*psi[n])(idx(_i,_j,_k)) - \
       (::min(real_t(0),C_adf_y(idx(_i,_j - grid->m_half,_k)))) * (*psi[n])(idx(_i,_j,_k)) + \
       (::max(real_t(0),C_adf_z(idx(_i,_j,_k + grid->p_half)))) * (*psi[n])(idx(_i,_j,_k)) - \
       (::min(real_t(0),C_adf_z(idx(_i,_j,_k - grid->m_half)))) * (*psi[n])(idx(_i,_j,_k))   \
     )
    C_mon_x(idx(iv + grid->p_half,j,k)) = C_adf_x(idx(iv + grid->p_half,j,k)) * where(
      C_adf_x(idx(iv + grid->p_half,j,k)) > 0,
      ::min(
        1, ::min(
          mpdata_beta_dn(iv  , j, k), 
          mpdata_beta_up(iv+1, j, k)
        )
      ),
      ::min(
        1, ::min(
          mpdata_beta_up(iv  , j, k), 
          mpdata_beta_dn(iv+1, j, k)
        )
      )
    );
  }
};
#endif
