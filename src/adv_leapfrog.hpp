/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @brief definition of the adv_leapfrog class with an implementation of the leapfrog advection scheme
 */
#ifndef ADV_LEAPFROG_HPP
#  define ADV_LEAPFROG_HPP

#  include "adv.hpp"
#  include "grd_arakawa-c-lorenz.hpp"

/// @brief an implementation of the leapfrog advection scheme
template <typename real_t> 
class adv_leapfrog : public adv<real_t> 
{
  public: const int stencil_extent() { return 3; }
  public: const int time_levels() { return 3; }
 
  public: const real_t courant_max() { return 1.; }
  public: const real_t courant_min() { return .5; } ; //TODO Wicker Skamarock 2002 
 
  private: grd_arakawa_c_lorenz<real_t> *grid;
  private: unique_ptr<adv_upstream<real_t>> fallback;

  public: adv_leapfrog(grd_arakawa_c_lorenz<real_t> *grid)
    : grid(grid), fallback(new adv_upstream<real_t>(grid))
  { 
    assert(this->stencil_extent() >= fallback->stencil_extent());
    assert(this->num_sclr_caches() == fallback->num_sclr_caches());
    assert(this->num_vctr_caches() == fallback->num_vctr_caches());
    assert(this->num_steps() == fallback->num_steps());
  }

  public: class op3D : public adv<real_t>::op3D
  {
    // private nested class for storing indices
    private:
    template <class idx>
    class indices
    {
      public: idx i_j_k, iph_j_k, imh_j_k, ip1_j_k, im1_j_k;
      public: indices(const mtx::idx &ijk, const grd_arakawa_c_lorenz<real_t> &grid) :
        i_j_k(  idx(ijk.i,               ijk.j, ijk.k)),
        iph_j_k(idx(ijk.i + grid.m_half, ijk.j, ijk.k)),
        imh_j_k(idx(ijk.i - grid.p_half, ijk.j, ijk.k)),
        ip1_j_k(idx(ijk.i + 1,           ijk.j, ijk.k)),
        im1_j_k(idx(ijk.i - 1,           ijk.j, ijk.k))
      { }
    };

    // private members 
    private: mtx::idx ijk;
    private: unique_ptr<typename adv_upstream<real_t>::op3D> flbkop;
    private: indices<mtx::idx_ijk> idxx;
    private: indices<mtx::idx_jki> idxy;
    private: indices<mtx::idx_kij> idxz;

    // ctor (initialisation of the indices)
    public: op3D(
      const mtx::idx &ijk, 
      const grd_arakawa_c_lorenz<real_t> &grid, 
      typename adv_upstream<real_t>::op3D *fallback
    ) : 
      ijk(ijk), idxx(ijk, grid), idxy(ijk, grid), idxz(ijk, grid)
    { 
      flbkop.reset(fallback);
    }
 
    // () operator 
    public: virtual void operator()(
      mtx::arr<real_t> *psi[], 
      const int n,  
      const int s,
      const mtx::arr<real_t> * const Cx, 
      const mtx::arr<real_t> * const Cy, 
      const mtx::arr<real_t> * const Cz
    ) 
    {
      *psi[n+1] = *psi[0]; // TODO: move from here to the solver?
      if (n == 1) // the leapfrog scheme
      {
        if (ijk.i_spans) op1D(psi, idxx, n, s, Cx, Cy, Cz);
        if (ijk.j_spans) op1D(psi, idxy, n, s, Cy, Cz, Cx);
        if (ijk.k_spans) op1D(psi, idxz, n, s, Cz, Cx, Cy); 
      }
      else if (n == 0) // upstream as a fallback scheme  
      {
        // fallback for the first timestep
        (*flbkop)(psi, n, s, Cx, Cy, Cz);
      }
      else assert(false);
    }

    // private methods
    private: 
    template <class indices>
    void op1D(
      mtx::arr<real_t>* psi[], 
      const indices &idx,
      const int n, 
      const int step,
      const mtx::arr<real_t> * const Cx, 
      const mtx::arr<real_t> * const, 
      const mtx::arr<real_t> * const 
    )
    {
      assert(step == 1);
      /// \f$ 
      ///   \psi^{n+1}_i = \psi^{n-1}_i - C^{n}_{i} \cdot (\psi^{n}_{i+1} - \psi^{n}_{i-1}) 
      /// \f$
      /// where C is the average Courant number for Arakawa C grid: 
      /// \f$ 
      ///   C^{n}_i=0.5\cdot(C^{n}_{i+1/2} + C^{n}_{i-1/2}) 
      /// \f$
      (*psi[n+1])(idx.i_j_k) -= 
        .5 * ( // average Courant number
          (*Cx)(idx.iph_j_k) + 
          (*Cx)(idx.imh_j_k)
        ) * (  
          (*psi[n])(idx.ip1_j_k) - 
          (*psi[n])(idx.im1_j_k) 
        );
    } 
  };

  public: virtual op3D *factory(
    const mtx::idx &ijk,
    mtx::arr<real_t> **, 
    mtx::arr<real_t> **
  ) 
  {
    return new op3D(ijk, *grid, fallback->factory(ijk, NULL, NULL));
  }
};
#endif
