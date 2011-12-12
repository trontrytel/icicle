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
  public: virtual const int num_vctr_caches() { return cache ? 3 : 0; }

  private: int iord;
  private: bool cache;
  private: grd_arakawa_c_lorenz<real_t> *grid;

  public: adv_mpdata(grd_arakawa_c_lorenz<real_t> *grid, int iord, bool cache = false) 
    : iord(iord), cache(cache), grid(grid)
  {
    if (iord <= 0) error_macro("iord (the number of iterations) must be > 0")
  }

    // using preprocessor macros as it's tricky make methods return parts of Blitz expressions 
#    define mpdata_frac(num, den) (where(den > real_t(0), (num) / (den), real_t(0)))
#    define mpdata_F(p1, p2, U) (real_t(.5) * (U + abs(U)) * p1 + real_t(.5) * (U - abs(U)) * p2)
#    define mpdata_A(pr, pl) mpdata_frac(pr - pl, pr + pl) 
#    define mpdata_B(pru, plu, prd, pld) (real_t(.5) * mpdata_frac(pru + plu - prd - pld, pru + plu + prd + pld))
#    define mpdata_V(Vru, Vlu, Vrd, Vld) (real_t(.25) * (Vru + Vlu + Vrd + Vld))
#    define mpdata_CA(pr, pl, U) ((abs(U) - pow(U,2)) * mpdata_A(pr, pl))
#    define mpdata_CB(pru, plu, prd, pld, U, V) (U * V * mpdata_B(pru, plu, prd, pld)) 

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
#    define mpdata_U(il,ic,ir) (\
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
        Range ir = iv+1, ic = iv + grid->p_half, il=iv;
        *tmp_v[dim] = mpdata_U(il, ic, ir);
        (*psi[n+1])(idx(i,j,k)) -= (
          mpdata_F((*psi[n])(idx(i  ,j,k)), (*psi[n])(idx(i+1,j,k)), (*tmp_v[dim])(idx(i + grid->p_half,j,k))) - 
          mpdata_F((*psi[n])(idx(i-1,j,k)), (*psi[n])(idx(i  ,j,k)), (*tmp_v[dim])(idx(i - grid->m_half,j,k)))
        );
      }
      else 
      {
        (*psi[n+1])(idx(i,j,k)) -= mpdata_F((*psi[n])(idx(i,  j,k)), (*psi[n])(idx(i+1,j,k)), mpdata_U(i  ,i+grid->p_half,i+1));
        (*psi[n+1])(idx(i,j,k)) += mpdata_F((*psi[n])(idx(i-1,j,k)), (*psi[n])(idx(i  ,j,k)), mpdata_U(i-1,i-grid->m_half,i  ));
        // TODO: check why the version below fails with a segfault in:
        //       blitz::FastArrayIteratorBase<float, 3, blitz::Array<float, 3> const&>::FastArrayIteratorBase ()
        //(*psi[n+1])(idx(i,j,k)) -= (
        //  mpdata_F((*psi[n])(idx(i,  j,k)), (*psi[n])(idx(i+1,j,k)), mpdata_U(i  ,i+grid->p_half,i+1)) -
        //  mpdata_F((*psi[n])(idx(i-1,j,k)), (*psi[n])(idx(i  ,j,k)), mpdata_U(i-1,i-grid->m_half,i  ))
        //);
      }
    }
#    undef mpdata_V
  }
#    undef mpdata_A
#    undef mpdata_B
#    undef mpdata_CA
#    undef mpdata_CB
#    undef mpdata_F
#  include "adv_hack.cpp"
};
#endif
