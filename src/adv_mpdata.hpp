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

#  include "adv.hpp"
#  include "grd_arakawa-c-lorenz.hpp"

template <typename real_t> 
class adv_mpdata : public adv<real_t> 
{
  public: const int stencil_extent() { return 3; }
  public: const int time_levels() { return 2; }
  public: const int num_steps() { return iord; }
  public: virtual const int num_vctr_caches() 
  { 
    return cache ? 3 : 0; 
  }

  private: int iord;
  private: bool cache;
  private: grd_arakawa_c_lorenz<real_t> *grid;

  public: adv_mpdata(grd_arakawa_c_lorenz<real_t> *grid, int iord, bool cache = false) 
    : iord(iord), cache(cache), grid(grid)
  {
    if (iord <= 0) error_macro("iord (the number of iterations) must be > 0")
  }

    // preprocessor macros are the only option as methods cannot return parts of Blitz expressions 
#    define mpdata_F(p1, p2, U) (.5 * (U + abs(U)) * p1 + .5 * (U - abs(U)) * p2)
#    define mpdata_A(pr, pl) \
       where(pr + pl > 0, \
         (pr - pl) / (pr + pl), \
         real_t(0.) \
       )
#    define mpdata_B(pru, plu, prd, pld) \
       where(pru + plu + prd + pld > 0, \
         .5 * (pru + plu - prd - pld) / (pru + plu + prd + pld), \
         real_t(0.) \
       )
#    define mpdata_V(Vru, Vlu, Vrd, Vld) (.25 * (Vru + Vlu + Vrd + Vld))
#    define mpdata_CA(pr, pl, U) ((abs(U) - pow(U,2)) * mpdata_A(pr, pl))
#    define mpdata_CB(pru, plu, prd, pld, U, V) (U * V * mpdata_B(pru, plu, prd, pld)) 
#    define mpdata_U (\
          mpdata_CA( \
            (*psi[n])(idx(ir, j, k)), (*psi[n])(idx(il, j, k)), /* pl, pr */ \
            Cx(idx(ic, j, k)) \
          ) - \
          mpdata_CB( \
            (*psi[n])(idx(ir, j+1, k)), (*psi[n])(idx(il, j+1, k)), /* pru, plu */ \
            (*psi[n])(idx(ir, j-1, k)), (*psi[n])(idx(il, j-1, k)), /* prd, pld */ \
            Cx(idx(ic, j, k)), mpdata_V( \
              Cy(idx(ir, j + grid->p_half, k)), Cy(idx(il, j + grid->p_half, k)), /* Vru, Vlu */ \
              Cy(idx(ir, j - grid->m_half, k)), Cy(idx(il, j - grid->m_half, k))  /* Vrd, Vld */ \
            ) \
          ) - \
          mpdata_CB( \
            (*psi[n])(idx(ir, j, k+1)), (*psi[n])(idx(il, j, k+1)), /* pru, plu */ \
            (*psi[n])(idx(ir, j, k-1)), (*psi[n])(idx(il, j, k-1)), /* prd, pld */ \
            Cx(idx(ic, j, k)), mpdata_V( \
              Cz(idx(ir, j, k + grid->p_half)), Cz(idx(il, j, k + grid->p_half)), /* Vru, Vlu */ \
              Cz(idx(ir, j, k - grid->m_half)), Cz(idx(il, j, k - grid->m_half))  /* Vrd, Vld */ \
            ) \
          )  \
        )
  public: 
  template <class idx>
  void op_helper(const real_t sign, const Range &il, const Range &ic, const Range &ir,
    const Range &i, const Range &j, const Range &k, 
    Array<real_t, 3> *psi[], const int n,
    const Array<real_t, 3> &Cx, 
    const Array<real_t, 3> &Cy, 
    const Array<real_t, 3> &Cz,
    Array<real_t, 3> *tmp[]
  )
  {
    (*psi[n+1])(idx(i,j,k)) += sign * (
      mpdata_F(
        (*psi[n])(idx(il, j, k)), (*psi[n])(idx(ir, j, k)), mpdata_U
      )
    );
  }

  public: 
  template <class idx>
  void op(int dim,
    Array<real_t, 3>* psi[], 
    Array<real_t, 3>* [], 
    Array<real_t, 3>* tmp_v[], 
    const Range &i, const Range &j, const Range &k, 
    const int n, const int step,
    const Array<real_t, 3> &Cx, 
    const Array<real_t, 3> &Cy, 
    const Array<real_t, 3> &Cz
  )
  {
    if (step == 1)
      (*psi[n+1])(idx(i,j,k)) -= (
        mpdata_F((*psi[n])(idx(i  ,j,k)), (*psi[n])(idx(i+1,j,k)), Cx(idx(i + grid->p_half,j,k))) - 
        mpdata_F((*psi[n])(idx(i-1,j,k)), (*psi[n])(idx(i  ,j,k)), Cx(idx(i - grid->m_half,j,k)))
      );
    else 
    {
      if (cache)
      {
        Range iv(i.first()-1, i.last());
        Range ir = iv+1, ic = iv+grid->p_half, il=iv;
        *tmp_v[dim] = mpdata_U;
        (*psi[n+1])(idx(i,j,k)) -= (
          mpdata_F((*psi[n])(idx(i  ,j,k)), (*psi[n])(idx(i+1,j,k)), (*tmp_v[dim])(idx(i + grid->p_half,j,k))) - 
          mpdata_F((*psi[n])(idx(i-1,j,k)), (*psi[n])(idx(i  ,j,k)), (*tmp_v[dim])(idx(i - grid->m_half,j,k)))
        );
      }
      else 
      {
        op_helper<idx>(-1, i  , i+grid->p_half, i+1, i, j, k, psi, n, Cx, Cy, Cz, tmp_v);
        op_helper<idx>(+1, i-1, i-grid->m_half, i  , i, j, k, psi, n, Cx, Cy, Cz, tmp_v);
      }
    }
  }
#    undef mpdata_A
#    undef mpdata_B
#    undef mpdata_V
#    undef mpdata_CA
#    undef mpdata_CB
#    undef mpdata_F
#  include "adv_hack.cpp"
};
#endif
