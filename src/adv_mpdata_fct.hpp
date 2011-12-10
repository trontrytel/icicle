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

  private: grd<real_t> *grid;
  public: adv_mpdata_fct(grd_arakawa_c_lorenz<real_t> *grid, int iord) 
    : adv_mpdata<real_t>(grid, iord, true), grid(grid)
  { }

/* Aniu! Tak jak obiecalem, nic nie zmieniam :)

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

#  define neighbours \
        (*psi[n+1])(i  ,j  ,k  ), \
        (*psi[n+1])(i-1,j  ,k  ), \
        (*psi[n+1])(i+1,j  ,k  ), \
        (*psi[n+1])(i  ,j-1,k  ), \
        (*psi[n+1])(i  ,j+1,k  ), \
        (*psi[n+1])(i  ,j  ,k-1), \
        (*psi[n+1])(i  ,j  ,k+1)
    if (step == 1)
    {
      psi_min(i,j,k) = min(neighbours); 
      psi_max(i,j,k) = max(neighbours); 
    }
    else
    {
      // calculating psi_min and psi_max
      psi_min(i,j,k) = min(psi_min(i,j,k), neighbours); 
      psi_max(i,j,k) = max(psi_max(i,j,k), neighbours); 

      // calculating Cx_mon, Cy_mon, Cz_mon
      fct_helper<idx_ijk>(psi, tmp_s, tmp_v, C_mon(0), C_adf(0), C_adf(1), C_adf(2), i, j, k, n);
      fct_helper<idx_jki>(psi, tmp_s, tmp_v, C_mon(1), C_adf(1), C_adf(2), C_adf(0), j, k, i, n);
      fct_helper<idx_kij>(psi, tmp_s, tmp_v, C_mon(2), C_adf(2), C_adf(0), C_adf(1), k, i, j, n);

      adv_mpdata<real_t>::op3D(psi, tmp_s, tmp_v, i, j, k, n, step, C_mon(0), C_mon(1), C_mon(2));
    }
#  undef neighbours

  }

  private:
  template <class idx>
  void fct_helper(
    Array<real_t, 3> *psi[], 
    Array<real_t, 3> *tmp_s[], 
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

#    define mpdata_frac(num, den) (where(den>real_t(0),(num)/(den),real_t(0)))  // TODO: move to mpdata
#    define mpdata_beta_up(i, j, k) mpdata_frac( \
       psi_max(idx(i,j,k)) - (*psi[n])(idx(i,j,k)), \
       (*psi[n+1])(max(0,C_adf_x(idx(i - grid->m_half,j,k)))) * (*psi[n])(idx(i-1,j,k)) - \
       (*psi[n+1])(min(0,C_adf_x(idx(i + grid->p_half,j,k)))) * (*psi[n])(idx(i+1,j,k)) + \
       (*psi[n+1])(max(0,C_adf_y(idx(i - grid->m_half,j,k)))) * (*psi[n])(idx(i-1,j,k)) - \
       (*psi[n+1])(min(0,C_adf_y(idx(i + grid->p_half,j,k)))) * (*psi[n])(idx(i+1,j,k)) + \
       (*psi[n+1])(max(0,C_adf_z(idx(i - grid->m_half,j,k)))) * (*psi[n])(idx(i-1,j,k)) - \
       (*psi[n+1])(min(0,C_adf_z(idx(i + grid->p_half,j,k)))) * (*psi[n])(idx(i+1,j,k))   \
     )
#    define mpdata_beta_dn(i, j, k) mpdata_frac(\
       (*psi[n])(idx(i,j,k)) - psi_min(idx(i,j,k)), \
       (*psi[n+1])(max(0,C_adf_x(idx(i + grid->p_half,j,k)))) * (*psi[n])(idx(i,j,k)) - \
       (*psi[n+1])(min(0,C_adf_x(idx(i - grid->m_half,j,k)))) * (*psi[n])(idx(i,j,k)) + \
       (*psi[n+1])(max(0,C_adf_y(idx(i + grid->p_half,j,k)))) * (*psi[n])(idx(i,j,k)) - \
       (*psi[n+1])(min(0,C_adf_y(idx(i - grid->m_half,j,k)))) * (*psi[n])(idx(i,j,k)) + \
       (*psi[n+1])(max(0,C_adf_z(idx(i + grid->p_half,j,k)))) * (*psi[n])(idx(i,j,k)) - \
       (*psi[n+1])(min(0,C_adf_z(idx(i - grid->m_half,j,k)))) * (*psi[n])(idx(i,j,k))   \
     )
    // 
    C_mon_x(idx(iv,j,k)) = C_adf_x(idx(iv,j,k)) * where(
      C_adf_x(idx(iv,j,k)) > 0,
      min(
        1, 
        mpdata_beta_dn(iv  , j, k), 
        mpdata_beta_up(iv+1, j, k)
      ),
      min(
        1, 
        mpdata_beta_up(iv  , j, k), 
        mpdata_beta_dn(iv+1, j, k)
      )
    );
    #undef mpdata_beta_up
    #undef mpdata_beta_dn
  }
*/
#  undef C_mon
#  undef C_adf
#  undef psi_min
#  undef psi_max
};
#endif
