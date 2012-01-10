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
  private: auto_ptr<rng> i, j, k;
  private: out<real_t> *output;
  private: int nx, ny, nz, halo;

  private: auto_ptr<arr<real_t> > Cx, Cy, Cz;

  private: auto_ptr<arr<real_t> > *psi_guard;
  private: arr<real_t> **psi; 

  private: auto_ptr<tmp<real_t> > cache; 

  public: slv_serial(stp<real_t> *setup,
    int i_min, int i_max, int nx,
    int j_min, int j_max, int ny,
    int k_min, int k_max, int nz, 
    quantity<si::time, real_t> dt // TODO: dt should not be needed here!
  )
    : fllbck(setup->fllbck), advsch(setup->advsch), output(setup->output), nx(nx), ny(ny), nz(nz)
  {
    if (fllbck != NULL) assert(advsch->stencil_extent() >= fllbck->stencil_extent());
    halo = (advsch->stencil_extent() - 1) / 2;

    // sanity checks
    {
      int len;
      if (halo > (len = i_max - i_min + 1) && nx != 1) 
        error_macro("halo length (" << halo << ") may not exceed domain extent in X (" << len << ")")
      if (halo > (len = j_max - j_min + 1) && ny != 1)
        error_macro("halo length (" << halo << ") may not exceed domain extent in Y (" << len << ")")
      if (halo > (len = k_max - k_min + 1) && nz != 1)
        error_macro("halo length (" << halo << ") may not exceed domain extent in Z (" << len << ")")
    }

    // memory allocation
    psi_guard = new auto_ptr<arr<real_t> >[advsch->time_levels()];
    psi = new arr<real_t>*[advsch->time_levels()];
    for (int n=0; n < advsch->time_levels(); ++n) 
    {
      psi_guard[n].reset(new arr<real_t>(
        setup->grid->rng_sclr(i_min, i_max, halo),
        setup->grid->rng_sclr(j_min, j_max, halo),
        setup->grid->rng_sclr(k_min, k_max, halo)
      ));
      psi[n] = psi_guard[n].get();
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
    i.reset(new rng(i_min, i_max));
    j.reset(new rng(j_min, j_max));
    k.reset(new rng(k_min, k_max));

    // initial condition
    setup->grid->populate_scalar_field(*i, *j, *k, psi[0], setup->intcond);

    // periodic boundary in all directions
    for (int s=this->first; s <= this->last; ++s) 
      this->hook_neighbour(s, this);
 
    // velocity fields
    {
      // TODO: this is grid-related logic! - to be moved from here
      rng 
        ix = setup->grid->rng_vctr(i->first(), i->last(), halo),
        jx = setup->grid->rng_sclr(j->first(), j->last(), halo),
        kx = setup->grid->rng_sclr(k->first(), k->last(), halo);
      Cx.reset(new arr<real_t>(ix, jx, kx));
      rng 
        iy = setup->grid->rng_sclr(i->first(), i->last(), halo),
        jy = setup->grid->rng_vctr(j->first(), j->last(), halo),
        ky = setup->grid->rng_sclr(k->first(), k->last(), halo);
      Cy.reset(new arr<real_t>(iy, jy, ky));
      rng 
        iz = setup->grid->rng_sclr(i->first(), i->last(), halo),
        jz = setup->grid->rng_sclr(j->first(), j->last(), halo),
        kz = setup->grid->rng_vctr(k->first(), k->last(), halo);
      Cz.reset(new arr<real_t>(iz, jz, kz));

      setup->grid->populate_courant_fields(
        ix, jx, kx, 
        iy, jy, ky, 
        iz, jz, kz, 
        Cx.get(), Cy.get(), Cz.get(), 
        setup->velocity, dt, 
        nx, ny, nz
      );
    }
  }

  public: ~slv_serial()
  {
    delete[] psi;
    delete[] psi_guard;
  }

  public: void record(const int n, const unsigned long t)
  {
    output->record(psi[n], *i, *j, *k, t);
  }

  public: typename arr<real_t>::arr_ret data(int n, const idx &idx)
  { 
    return (*psi[n])(idx);
  }

  public: void fill_halos(int n)
  {
    fill_halos_helper<idx_ijk>(psi, slv<real_t>::left, n, i->first() - halo, i->first() - 1, *j, *k, nx);
    fill_halos_helper<idx_ijk>(psi, slv<real_t>::rght, n, i->last() + 1, i->last() + halo, *j, *k, nx);
    fill_halos_helper<idx_jki>(psi, slv<real_t>::fore, n, j->first() - halo, j->first() - 1, *k, *i, ny);
    fill_halos_helper<idx_jki>(psi, slv<real_t>::hind, n, j->last() + 1, j->last() + halo, *k, *i, ny);
    fill_halos_helper<idx_kij>(psi, slv<real_t>::base, n, k->first() - halo, k->first() - 1, *i, *j, nz);
    fill_halos_helper<idx_kij>(psi, slv<real_t>::apex, n, k->last() + 1, k->last() + halo, *i, *j, nz);
  }

  private:
  template<class idx>
  void fill_halos_helper(arr<real_t> *psi[], int nghbr, int n, 
    int i_min, int i_max, const rng &j, const rng &k, int mod
  )
  {
    if (mod == 1)
      for (int ii = i_min; ii <= i_max; ++ii)
        (*psi[n])(idx(rng(ii), j, k)) =
          this->nghbr_data(nghbr, n, idx(rng(0,0), j, k)); // only happens with periodic boundary
    else 
      (*psi[n])(idx(rng(i_min, i_max), j, k)) =
        this->nghbr_data(nghbr, n, idx(rng((i_min + mod) % mod, (i_max + mod) % mod), j, k));
  }

  public: void advect(adv<real_t> *a, int n, int s, quantity<si::time, real_t> dt)
  {
    a->op3D(psi, cache->sclr, cache->vctr, *i, *j, *k, n, s, Cx.get(), Cy.get(), Cz.get());
  }

  public: void cycle_arrays(const int n)
  {
    switch (advsch->time_levels())
    {
      case 2: // e.g. upstream, mpdata
        cycleArrays(*psi[n], *psi[n+1]); 
        break;
      case 3: // e.g. leapfrog
        cycleArrays(*psi[n-1], *psi[n], *psi[n+1]); 
        break;
      default: assert(false);
    }
  }

  public: void integ_loop(unsigned long nt, quantity<si::time, real_t> dt) 
  {
    assert(false); // TODO: this should not be needed with proper class hierarchy
  }
};
#endif
