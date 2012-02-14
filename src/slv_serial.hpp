/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011 - February 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef SLV_SERIAL_HPP
#  define SLV_SERIAL_HPP

#  include "slv.hpp"
#  include "stp.hpp"
#  include "tmp.hpp"

template <typename real_t>
class slv_serial : public slv<real_t>
{
  private: adv<real_t> *fllbck, *advsch;
  private: auto_ptr<mtx::idx> ijk;
  private: out<real_t> *output;
  private: stp<real_t> *setup;
  private: vector<int> halo_sclr;
  private: int halo_vctr;

  private: auto_ptr<mtx::arr<real_t> > Cx, Cy, Cz;
  private: mtx::arr<real_t> **Qx, **Qy, **Qz;
  private: vector<ptr_vector<mtx::arr<real_t> > > psi;
  private: auto_ptr<tmp<real_t> > cache; 

  public: slv_serial(stp<real_t> *setup, out<real_t> *output,
    int i_min, int i_max,
    int j_min, int j_max,
    int k_min, int k_max
  )
    : fllbck(setup->fllbck), advsch(setup->advsch), output(output), setup(setup)
  {
    // computing required halo lengths
    if (fllbck != NULL) assert(advsch->stencil_extent() >= fllbck->stencil_extent());
    halo_vctr = (advsch->stencil_extent() - 1) / 2;
    halo_sclr.assign(setup->equations->n_vars(), halo_vctr);
    if (!setup->velocity->is_constant())
    {
      // enlarged halo needed for t=n+1/2 velocity extrapolation
      // TODO: stencil_extent()-like method in velocity?
      if (setup->equations->has_qx()) halo_sclr[setup->equations->idx_qx()] += 1; 
      if (setup->equations->has_qy()) halo_sclr[setup->equations->idx_qy()] += 1; 
      if (setup->equations->has_qz()) halo_sclr[setup->equations->idx_qz()] += 1; 
    }

    // memory allocation
    psi.resize(setup->equations->n_vars());
    for (int e = 0; e < setup->equations->n_vars(); ++e)
    {
      // sanity checks
      {
        int len;
        if (halo_sclr[e] > (len = i_max - i_min + 1) && setup->grid->nx() != 1) 
          error_macro("halo length (" << halo_sclr[e] << ") may not exceed domain extent in X (" << len << ")")
        if (halo_sclr[e] > (len = j_max - j_min + 1) && setup->grid->ny() != 1)
          error_macro("halo length (" << halo_sclr[e] << ") may not exceed domain extent in Y (" << len << ")")
        if (halo_sclr[e] > (len = k_max - k_min + 1) && setup->grid->nz() != 1)
          error_macro("halo length (" << halo_sclr[e] << ") may not exceed domain extent in Z (" << len << ")")
      }

      int tlevs = advsch->time_levels();
      if (fllbck != NULL) assert(tlevs > fllbck->time_levels()); // TODO: OK ?
      for (int n = 0; n < advsch->time_levels(); ++n) 
      {
        psi[e].push_back(new mtx::arr<real_t>(
          setup->grid->rng_sclr(i_min, i_max, j_min, j_max, k_min, k_max, halo_sclr[e])
        ));
      }
    }

    // momentum pointers
    if (setup->velocity->is_constant()) Qx = Qy = Qz = NULL; // TODO: czy to potrzebne?
    else 
    {
      if (setup->equations->has_qx()) Qx = psi[setup->equations->idx_qx()].c_array(); 
      if (setup->equations->has_qy()) Qy = psi[setup->equations->idx_qy()].c_array(); 
      if (setup->equations->has_qz()) Qz = psi[setup->equations->idx_qz()].c_array(); 
    }

    // caches
    assert(fllbck == NULL || fllbck->num_sclr_caches() == 0);
    assert(fllbck == NULL || fllbck->num_vctr_caches() == 0);
    cache.reset(new tmp<real_t>(advsch->num_vctr_caches(), advsch->num_sclr_caches(), setup->grid, 
      halo_vctr, // halo_sclr > halo_vctr only for non-const velocities, the cache is only ised in adv_*
      i_min, i_max,
      j_min, j_max,
      k_min, k_max
    ));

    // indices
    ijk.reset(new mtx::idx_ijk(
      mtx::rng(i_min, i_max), 
      mtx::rng(j_min, j_max), 
      mtx::rng(k_min, k_max)
    ));

    // initial condition
    for (int e = 0; e < setup->equations->n_vars(); ++e)
      setup->intcond->populate_scalar_field(setup->equations->var_name(e), *ijk, psi[e][0]); 

    // periodic boundary in all directions
    for (int s=this->first; s <= this->last; ++s) 
      this->hook_neighbour(s, this);
 
    // velocity fields
    Cx.reset(new mtx::arr<real_t>(setup->grid->rng_vctr_x(*ijk, halo_vctr)));
    Cy.reset(new mtx::arr<real_t>(setup->grid->rng_vctr_y(*ijk, halo_vctr)));
    Cz.reset(new mtx::arr<real_t>(setup->grid->rng_vctr_z(*ijk, halo_vctr)));
    if (setup->velocity->is_constant()) 
    {
      // filling Cx, Cy, Cz with constant velocity field
      setup->velocity->populate_courant_fields(-1, // TODO: some nicer calling sequence?
        Cx.get(), Cy.get(), Cz.get(), setup->dt, NULL, NULL, NULL // ... Qx, Qy, Qx should anyhow contain NULLs
      );
    }
    else
    {
      // first guess for velocity fields assuming constant momenta
      // TODO two or more, use a for!
      int n = 0;
      if (setup->equations->has_qx()) 
      {
        fill_halos(setup->equations->idx_qx(), n); 
        *Qx[n+1] = *Qx[n];
      }
      if (setup->equations->has_qy()) 
      {
        fill_halos(setup->equations->idx_qy(), n); 
        *Qy[n+1] = *Qy[n];
      }
      if (setup->equations->has_qz()) 
      {
        fill_halos(setup->equations->idx_qz(), n); 
        *Qz[n+1] = *Qz[n]; 
      }
      setup->velocity->populate_courant_fields(n + 1,
        Cx.get(), Cy.get(), Cz.get(), setup->dt, Qx, Qy, Qz
      );
#  ifndef NDEBUG
      if (setup->equations->has_qx()) Qx[n+1]->fill_with_nans();
      if (setup->equations->has_qy()) Qy[n+1]->fill_with_nans();
      if (setup->equations->has_qz()) Qz[n+1]->fill_with_nans();
#  endif
    }
  }

  public: void record(const int n, const unsigned long t)
  {
    for (int e = 0; e < setup->equations->n_vars(); ++e)
      output->record(e, psi[e][n], *ijk, t);
  }

  public: typename mtx::arr<real_t>::type data(int e, int n, const mtx::idx &idx)
  { 
    return psi[e][n](idx);
  }

  public: void fill_halos(const int e, const int n)
  {
    mtx::rng 
      i_all(ijk->lbound(mtx::i) - halo_sclr[e], ijk->ubound(mtx::i) + halo_sclr[e]), // TODO: replace with psi[e][n]->i ?
      j_all(ijk->lbound(mtx::j) - halo_sclr[e], ijk->ubound(mtx::j) + halo_sclr[e]); // TODO: replace with psi[e][n]->j ?
    fill_halos_helper<mtx::idx_ijk>(slv<real_t>::left, e, n, ijk->lbound(mtx::i) - halo_sclr[e], ijk->lbound(mtx::i) - 1,            ijk->j, ijk->k, setup->grid->nx());
    fill_halos_helper<mtx::idx_ijk>(slv<real_t>::rght, e, n, ijk->ubound(mtx::i) + 1,            ijk->ubound(mtx::i) + halo_sclr[e], ijk->j, ijk->k, setup->grid->nx());
    fill_halos_helper<mtx::idx_jki>(slv<real_t>::fore, e, n, ijk->lbound(mtx::j) - halo_sclr[e], ijk->lbound(mtx::j) - 1,            ijk->k, i_all,  setup->grid->ny());
    fill_halos_helper<mtx::idx_jki>(slv<real_t>::hind, e, n, ijk->ubound(mtx::j) + 1,            ijk->ubound(mtx::j) + halo_sclr[e], ijk->k, i_all,  setup->grid->ny());
    fill_halos_helper<mtx::idx_kij>(slv<real_t>::base, e, n, ijk->lbound(mtx::k) - halo_sclr[e], ijk->lbound(mtx::k) - 1,            i_all,  j_all,  setup->grid->nz());
    fill_halos_helper<mtx::idx_kij>(slv<real_t>::apex, e, n, ijk->ubound(mtx::k) + 1,            ijk->ubound(mtx::k) + halo_sclr[e], i_all,  j_all,  setup->grid->nz());
  }

  private:
  template<class idx>
  void fill_halos_helper(const int nghbr, const int e, const int n, 
    const int i_min, const int i_max, const mtx::rng &j, const mtx::rng &k, int mod
  )
  {
    if (mod == 1)
    {
      for (int ii = i_min; ii <= i_max; ++ii)
        psi[e][n](idx(mtx::rng(ii), j, k)) =
          this->nghbr_data(nghbr, e, n, idx(mtx::rng(0,0), j, k)); // only happens with periodic boundary
    }
    else 
    {
      psi[e][n](idx(mtx::rng(i_min, i_max), j, k)) =
        this->nghbr_data(nghbr, e, n, idx(mtx::rng((i_min + mod) % mod, (i_max + mod) % mod), j, k));
    }
  }

  public: void advect(int e, int n, int s, adv<real_t> *a)
  {
    a->op3D(psi[e].c_array(), cache->sclr, cache->vctr, *ijk, n, s, 
      Cx.get(), Cy.get(), Cz.get() 
    );
  }

  public: void update_velocities(const int n) // TODO: rename to update courants?
  {
    assert(!setup->velocity->is_constant());
    assert(!setup->equations->has_qx() || Qx != NULL);
    assert(!setup->equations->has_qy() || Qy != NULL); 
    assert(!setup->equations->has_qz() || Qz != NULL); 
    setup->velocity->populate_courant_fields(n,
      Cx.get(), Cy.get(), Cz.get(), setup->dt, Qx, Qy, Qz
    );
  }

  public: void cycle_arrays(const int e, const int n)
  {
    switch (advsch->time_levels())
    {
      case 2: // e.g. upstream, mpdata
        mtx::cycle(psi[e][n], psi[e][n+1]);
        break;
      case 3: // e.g. leapfrog
        mtx::cycle(psi[e][n-1], psi[e][n], psi[e][n+1]);
        break;
      default: assert(false);
    }
  }

  public: void stash_cycle(int e, int n)
  {
    static mtx::arr<real_t> stash(psi[e][n]);
    // TODO: static size + assert() if the same size
    mtx::cycle(stash, psi[e][n]);
  }

  public: void integ_loop() 
  {
    assert(false); // TODO: this should not be needed with proper class hierarchy
  }
};
#endif
