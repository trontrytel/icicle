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
  adv_hack_macro // workaround for virtual template methods
  public: const int num_steps() { return iord; }
  public: virtual const int num_vctr_caches() { return 1; }

  private: int iord;
  private: grd_arakawa_c_lorenz<real_t> *grid;

  public: adv_mpdata(grd_arakawa_c_lorenz<real_t> *grid, int iord) 
    : iord(iord), grid(grid), adv_upstream<real_t>(grid)
  {
    if (iord <= 0) error_macro("iord (the number of iterations) must be > 0")
  }

    // using preprocessor macros as it's tricky make methods return parts of Blitz expressions 
#    define mpdata_frac(num, den) (where(den > real_t(0), (num) / (den), real_t(0)))
#    define mpdata_A(pr, pl) mpdata_frac(pr - pl, pr + pl) 
#    define mpdata_B(pru, plu, prd, pld) (real_t(.5) * mpdata_frac(pru + plu - prd - pld, pru + plu + prd + pld))
#    define mpdata_V(Vru, Vlu, Vrd, Vld) (real_t(.25) * (Vru + Vlu + Vrd + Vld))
#    define mpdata_CA(pr, pl, U) ((abs(U) - pow(U,2)) * mpdata_A(pr, pl))
#    define mpdata_CB(pru, plu, prd, pld, U, V) (U * V * mpdata_B(pru, plu, prd, pld)) 

  protected:
  template <class idx>
  void mpdata_U(
    arr<real_t> *C_adf,
    arr<real_t> *psi[], const int n,
    const Range &i, const Range &j, const Range &k,
    const arr<real_t> &Cx, const arr<real_t> &Cy, const arr<real_t> &Cz
  )
  {
    // instead of computing u_{i+1/2} and u_{i-1/2} for all i
    // we compute u_{i+1/2} for im=(i-1, ... i)
    Range im(i.first() - 1, i.last());
    Range ir = im + 1, ic = im + grid->p_half, il = im;
    (*C_adf)(idx(Range(i.first() - grid->m_half, i.last() + grid->p_half), j, k)) = (
      mpdata_CA( 
        (*psi[n])(idx(ir, j, k)), (*psi[n])(idx(il, j, k)), /* pl, pr */ 
        Cx(idx(ic, j, k)) 
      ) - 
      mpdata_CB( 
        (*psi[n])(idx(ir, j+1, k)), (*psi[n])(idx(il, j+1, k)), /* pru, plu */ 
        (*psi[n])(idx(ir, j-1, k)), (*psi[n])(idx(il, j-1, k)), /* prd, pld */ 
        Cx(idx(ic, j, k)), mpdata_V( 
          Cy(idx(ir, j + grid->p_half, k)), Cy(idx(il, j + grid->p_half, k)), /* Vru, Vlu */ 
          Cy(idx(ir, j - grid->m_half, k)), Cy(idx(il, j - grid->m_half, k))  /* Vrd, Vld */ 
        ) 
      ) - 
      mpdata_CB( 
        (*psi[n])(idx(ir, j, k+1)), (*psi[n])(idx(il, j, k+1)), /* pru, plu */ 
        (*psi[n])(idx(ir, j, k-1)), (*psi[n])(idx(il, j, k-1)), /* prd, pld */ 
        Cx(idx(ic, j, k)), mpdata_V( 
          Cz(idx(ir, j, k + grid->p_half)), Cz(idx(il, j, k + grid->p_half)), /* Vru, Vlu */ 
          Cz(idx(ir, j, k - grid->m_half)), Cz(idx(il, j, k - grid->m_half))  /* Vrd, Vld */ 
        ) 
      )  
    );
  }

  public: 
  template <class idx>
  void op(int dim,
    arr<real_t>* psi[], 
    arr<real_t>* [], 
    arr<real_t>* tmp_v[], 
    const Range &i, const Range &j, const Range &k, 
    const int n, const int step,
    const arr<real_t> &Cx, const arr<real_t> &Cy, const arr<real_t> &Cz
  )
  {
    if (step == 1) 
      adv_upstream<real_t>::template op<idx>(dim, psi, NULL, NULL, i, j, k, n, 1, Cx, Cy, Cz);
    else 
    {
      mpdata_U<idx>(tmp_v[0], psi, n, i, j, k, Cx, Cy, Cz);
      adv_upstream<real_t>::template op<idx>(dim, psi, NULL, NULL, i, j, k, n, 1, *tmp_v[0], *tmp_v[1], *tmp_v[2]);
      //                                                                                     ^^^^^^^^^^^^^^^^^^^^ TODO: these are non-existant!
    }
  }
};
#endif
