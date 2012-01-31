/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef DOM_SERIAL_HPP
#  define DOM_SERIAL_HPP

#  include "slv.hpp"
#  include "stp.hpp"
#  include "tmp.hpp"

template <typename real_t>
class slv_serial : public slv<real_t>
{
  private: adv<real_t> *fllbck, *advsch;
  private: auto_ptr<mtx::idx> ijk;
  private: auto_ptr<mtx::idx> ijk_all;
  private: out<real_t> *output;
  private: stp<real_t> *setup;
  private: int halo;

  private: auto_ptr<mtx::arr<real_t> > Cx, Cy, Cz;
  private: vector<ptr_vector<mtx::arr<real_t> > > psi;
  private: auto_ptr<tmp<real_t> > cache; 

  public: slv_serial(stp<real_t> *setup, out<real_t> *output,
    int i_min, int i_max,
    int j_min, int j_max,
    int k_min, int k_max
  )
    : fllbck(setup->fllbck), advsch(setup->advsch), output(output), setup(setup)
  {
    if (fllbck != NULL) assert(advsch->stencil_extent() >= fllbck->stencil_extent());
    halo = (advsch->stencil_extent() - 1) / 2;

    // sanity checks
    {
      int len;
      if (halo > (len = i_max - i_min + 1) && setup->grid->nx() != 1) 
        error_macro("halo length (" << halo << ") may not exceed domain extent in X (" << len << ")")
      if (halo > (len = j_max - j_min + 1) && setup->grid->ny() != 1)
        error_macro("halo length (" << halo << ") may not exceed domain extent in Y (" << len << ")")
      if (halo > (len = k_max - k_min + 1) && setup->grid->nz() != 1)
        error_macro("halo length (" << halo << ") may not exceed domain extent in Z (" << len << ")")
    }

    // memory allocation
    psi.resize(setup->equations->n_vars());
    for (int e = 0; e < setup->equations->n_vars(); ++e)
    {
      for (int n = 0; n < advsch->time_levels(); ++n) 
      {
        psi[e].push_back(new mtx::arr<real_t>(
          setup->grid->rng_sclr(i_min, i_max, j_min, j_max, k_min, k_max, halo)
        ));
      }
    }

    // caches
    assert(fllbck == NULL || fllbck->num_sclr_caches() == 0);
    assert(fllbck == NULL || fllbck->num_vctr_caches() == 0);
    cache.reset(new tmp<real_t>(advsch->num_vctr_caches(), advsch->num_sclr_caches(), setup->grid, halo,
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
    ijk_all.reset(new mtx::idx_ijk(
      mtx::rng(i_min - halo, i_max + halo),
      mtx::rng(j_min - halo, j_max + halo),
      mtx::rng(k_min - halo, k_max + halo)
    ));

    // initial condition
    for (int e = 0; e < setup->equations->n_vars(); ++e)
      setup->intcond->populate_scalar_field(setup->equations->var_name(e), *ijk, psi[e][0]); 

    // periodic boundary in all directions
    for (int s=this->first; s <= this->last; ++s) 
      this->hook_neighbour(s, this);
 
    // velocity fields
    Cx.reset(new mtx::arr<real_t>(setup->grid->rng_vctr_x(*ijk, halo)));
    Cy.reset(new mtx::arr<real_t>(setup->grid->rng_vctr_y(*ijk, halo)));
    Cz.reset(new mtx::arr<real_t>(setup->grid->rng_vctr_z(*ijk, halo)));
    setup->grid->populate_courant_fields(Cx.get(), Cy.get(), Cz.get(), setup->velocity, setup->dt);
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
    fill_halos_helper<mtx::idx_ijk>(slv<real_t>::left, e, n, ijk->lbound(mtx::i) - halo, ijk->lbound(mtx::i) - 1,    ijk->j, ijk->k, setup->grid->nx());
    fill_halos_helper<mtx::idx_ijk>(slv<real_t>::rght, e, n, ijk->ubound(mtx::i) + 1,    ijk->ubound(mtx::i) + halo, ijk->j, ijk->k, setup->grid->nx());
    fill_halos_helper<mtx::idx_jki>(slv<real_t>::fore, e, n, ijk->lbound(mtx::j) - halo, ijk->lbound(mtx::j) - 1,    ijk->k, ijk_all->i, setup->grid->ny());
    fill_halos_helper<mtx::idx_jki>(slv<real_t>::hind, e, n, ijk->ubound(mtx::j) + 1,    ijk->ubound(mtx::j) + halo, ijk->k, ijk_all->i, setup->grid->ny());
    fill_halos_helper<mtx::idx_kij>(slv<real_t>::base, e, n, ijk->lbound(mtx::k) - halo, ijk->lbound(mtx::k) - 1,    ijk_all->i, ijk_all->j, setup->grid->nz());
    fill_halos_helper<mtx::idx_kij>(slv<real_t>::apex, e, n, ijk->ubound(mtx::k) + 1,    ijk->ubound(mtx::k) + halo, ijk_all->i, ijk_all->j, setup->grid->nz());
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
    a->op3D(psi[e].c_array(), cache->sclr, cache->vctr, *ijk, n, s, Cx.get(), Cy.get(), Cz.get());
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

  public: void integ_loop() 
  {
    assert(false); // TODO: this should not be needed with proper class hierarchy
  }
};
#endif
