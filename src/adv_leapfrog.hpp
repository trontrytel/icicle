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
#  include "grd.hpp"

/// @brief an implementation of the leapfrog advection scheme
template <typename real_t> 
class adv_leapfrog : public adv<real_t> 
{
  public: const int stencil_extent() const { return 3; }
  public: const int time_levels() const { return 3; }
 
  public: const real_t courant_max() const { return 1.; }
  public: const real_t courant_min() const { return .5; } ; //TODO Wicker Skamarock 2002 
 
  private: unique_ptr<adv_upstream<real_t>> fallback;

  public: adv_leapfrog()
    : fallback(new adv_upstream<real_t>())
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
      public: const idx i_j_k, iph_j_k, imh_j_k, ip1_j_k, im1_j_k;
      public: indices(
        const mtx::rng &i, 
        const mtx::rng &j, 
        const mtx::rng &k
      ) :
        i_j_k(  idx(i,               j, k)),
        iph_j_k(idx(i + static_rational<1,2>(), j, k)),
        imh_j_k(idx(i - static_rational<1,2>(), j, k)),
        ip1_j_k(idx(i + 1,           j, k)),
        im1_j_k(idx(i - 1,           j, k))
      { }
    };

    // private members 
    private: unique_ptr<typename adv<real_t>::op3D> flbkop;
    private: const indices<mtx::idx_ijk> idxx;
    private: const indices<mtx::idx_jki> idxy;
    private: const indices<mtx::idx_kij> idxz;

    // ctor (initialisation of the indices)
    public: op3D(
      const mtx::idx &ijk, 
      typename adv<real_t>::op3D *fallback
    ) : 
      adv<real_t>::op3D(ijk),
      idxx(ijk.i, ijk.j, ijk.k), 
      idxy(ijk.j, ijk.k, ijk.i), 
      idxz(ijk.k, ijk.i, ijk.j)
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
        if (this->do_x()) op1D(psi, idxx, n, s, Cx, Cy, Cz);
        if (this->do_y()) op1D(psi, idxy, n, s, Cy, Cz, Cx);
        if (this->do_z()) op1D(psi, idxz, n, s, Cz, Cx, Cy); 
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

  public: typename adv<real_t>::op3D *factory(
    const mtx::idx &ijk,
    mtx::arr<real_t> **, 
    mtx::arr<real_t> **,
    bool 
  ) const
  {
    return new op3D(ijk, fallback->factory(ijk, NULL, NULL, false));
  }
};
#endif
