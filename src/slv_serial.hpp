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
  private: unique_ptr<mtx::idx> ijk;
  private: out<real_t> *output;
  private: stp<real_t> *setup;
  private: vector<int> halo_sclr;
  private: int halo_vctr;

  private: ptr_vector<mtx::arr<real_t>> C; // TODO: nullable? (uwaga na kod w FCT! ii = i+/-1
  private: vector<ptr_vector<mtx::arr<real_t>>> psi;
  private: unique_ptr<tmp<real_t>> cache; 
  private: ptr_vector<nullable<mtx::arr<real_t>>> rhs_R, rhs_C;

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
    halo_sclr.assign(setup->eqsys->n_vars(), halo_vctr);
    if (!setup->velocity->is_constant())
    {
      // enlarged halo needed for t=n+1/2 velocity extrapolation
      // TODO: stencil_extent()-like method in velocity?
      for (int e = 0; e < setup->eqsys->n_vars(); ++e)
        if (setup->eqsys->var_dynamic(e)) halo_sclr[e] += 1;
    }

    // memory allocation
    psi.resize(setup->eqsys->n_vars());
    for (int e = 0; e < setup->eqsys->n_vars(); ++e)
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

      if (!setup->eqsys->is_homogeneous(e))
      {
        rhs_R.push_back(new mtx::arr<real_t>(setup->grid->rng_sclr(
          i_min, i_max, j_min, j_max, k_min, k_max, 0)
        ));
      } else rhs_R.push_back(NULL);
    }

    // caches (TODO: move to adv.hpp...)
    assert(fllbck == NULL || fllbck->num_sclr_caches() == 0);
    assert(fllbck == NULL || fllbck->num_vctr_caches() == 0);
    cache.reset(new tmp<real_t>(advsch->num_vctr_caches(), advsch->num_sclr_caches(), setup->grid, 
      halo_vctr + (setup->velocity->is_constant() ? 0 : 1), /* TODO TODO TODO it should be max(halo_sclr, halo_vctr) */ // halo_sclr > halo_vctr only for non-const velocities, the cache is only used in adv_*
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
    {
      int n = 0;
      for (int e = 0; e < setup->eqsys->n_vars(); ++e)
        setup->intcond->populate_scalar_field(setup->eqsys->var_name(e), *ijk, psi[e][n]); 
    }

    // periodic boundary in all directions
    for (int s=this->first; s <= this->last; ++s) 
      this->hook_neighbour(s, this);
 
    // velocity fields 
    C.push_back(new mtx::arr<real_t>(setup->grid->rng_vctr_x(*ijk, halo_vctr)));
    C.push_back(new mtx::arr<real_t>(setup->grid->rng_vctr_y(*ijk, halo_vctr)));
    C.push_back(new mtx::arr<real_t>(setup->grid->rng_vctr_z(*ijk, halo_vctr)));
    if (setup->velocity->is_constant()) 
    {
      // filling C with constant velocity field
      setup->velocity->populate_courant_fields(-1,-1, // TODO: some nicer calling sequence?
        &C[0], &C[1], &C[2], setup->dt, NULL, NULL, NULL 
      );
    }
    //else 
    //{
    //  // filling C with zeros
    //  for (int xyz = 0; xyz < 3; ++xyz)
    //    C[xyz](C[xyz].ijk) = real_t(0);
    //}
  }

  public: void copy(int from, int to)
  {
    for (int e = 0; e < setup->eqsys->n_vars(); ++e)
      psi[e][to] = psi[e][from];
  }

  public: void fill_with_nans(int e, int n)
  { 
#  ifndef NDEBUG
    psi[e][n].fill_with_nans();
#  endif
  }

  public: void record(const int n, const unsigned long t)
  {
    for (int e = 0; e < setup->eqsys->n_vars(); ++e)
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
      &C[0], &C[1], &C[2] 
    );
  }

  public: void update_courants(const int g, const int nm1, const int nm0) // TODO: this method deserves a major rewrite!
  {
    assert(!setup->velocity->is_constant());

    mtx::arr<real_t> **Q[] = {NULL, NULL, NULL};
    for (int xyz = 0; xyz < 3; ++xyz)
    {
      if (setup->eqsys->velmap(g, xyz).size() > 0) // if we have any source terms
      { 
        Q[xyz] = psi[setup->eqsys->velmap(g, xyz).begin()->first].c_array(); 
        // TODO: something more universal would be better
        assert(setup->eqsys->velmap(g, xyz).begin()->second == 1); 
      }
    }

    // 
    for (int xyz = 0; xyz < 3; ++xyz) 
    {
      if (setup->eqsys->velmap(g, xyz).size() > 1) 
      {
        assert(Q[xyz] != NULL);
        for (int i = 1; i < setup->eqsys->velmap(g, xyz).size(); ++i) 
        {
          int 
            var = setup->eqsys->velmap(g, xyz)[i].first, 
            pow = setup->eqsys->velmap(g, xyz)[i].second;
          if (pow == -1) 
          {
            {
              assert(isfinite(sum((psi[var][nm0])(psi[var][nm0].ijk))));
              assert(isfinite(sum((psi[var][nm1])(psi[var][nm1].ijk))));

              (*(Q[xyz][nm0]))(Q[xyz][nm0]->ijk) = where(
                (psi[var][nm0])(psi[var][nm0].ijk) != 0,
                (*(Q[xyz][nm0]))(Q[xyz][nm0]->ijk) / (psi[var][nm0])(psi[var][nm0].ijk),
                real_t(0)
              ); // TODO
              (*(Q[xyz][nm1]))(Q[xyz][nm1]->ijk) = where(
                (psi[var][nm1])(psi[var][nm1].ijk) != 0,
                (*(Q[xyz][nm1]))(Q[xyz][nm1]->ijk) / (psi[var][nm1])(psi[var][nm1].ijk),
                real_t(0)
              ); // TODO
              //(*(Q[xyz][nm1]))(Q[xyz][nm1]->ijk) /= (psi[var][nm1])(psi[var][nm1].ijk); // TODO: two or more, use a for
            }
          }
          else assert(false);
        }
      }
    }
 
    // compute the velocities
    setup->velocity->populate_courant_fields(nm0, nm1,
      &C[0], &C[1], &C[2], setup->dt, Q[0], Q[1], Q[2]
    );

    // 
    for (int xyz = 0; xyz < 3; ++xyz) 
    {
      if (setup->eqsys->velmap(g, xyz).size() > 1) 
      {
        for (int i = 1; i < setup->eqsys->velmap(g, xyz).size(); ++i) 
        {
          int 
            var = setup->eqsys->velmap(g, xyz)[i].first, 
            pow = setup->eqsys->velmap(g, xyz)[i].second;
          if (pow == -1) 
          {
            {
              assert(isfinite(sum((psi[var][nm0])(psi[var][nm0].ijk))));
              (*(Q[xyz][nm0]))(Q[xyz][nm0]->ijk) *= (psi[var][nm0])(psi[var][nm0].ijk); // TODO
              assert(isfinite(sum((psi[var][nm1])(psi[var][nm1].ijk))));
              (*(Q[xyz][nm1]))(Q[xyz][nm1]->ijk) *= (psi[var][nm1])(psi[var][nm1].ijk); // TODO: two or more, use a for
            }
          }
          else assert(false);
        }
      }
    }

    // TODO: sanity checks for Courant limits
    cerr << "Courant x: ";
    cerr << min(C[0]) << " ... ";
    cerr << max(C[0]) << endl;
    //cerr << "Courant y: " << min(C[1]) << " ... " << max(C[1]) << endl;
    //cerr << "Courant z: " << min(C[2]) << " ... " << max(C[2]) << endl;
  }

  /// \param n the time level to use for updating the forcings
  private: mtx::arr<real_t> **tmpvec = NULL;
  public: void update_forcings(int n)
  {
    if (tmpvec == NULL) tmpvec = new mtx::arr<real_t>*[setup->eqsys->n_vars()]; // TODO: unique_ptr

    for (int e = 0; e < setup->eqsys->n_vars(); ++e) tmpvec[e] = psi[e].c_array()[n];

    for (int e = 0; e < setup->eqsys->n_vars(); ++e)
    {
      if (!setup->eqsys->is_homogeneous(e))
      {
        (rhs_R[e])(rhs_R[e].ijk) = real_t(0);
        for (int i = 0; i < setup->eqsys->var_n_rhs_terms(e); ++i)
        {
           setup->eqsys->var_rhs_term(e, i).explicit_part(rhs_R[e], tmpvec, *ijk);
           assert(isfinite(sum((rhs_R[e])(rhs_R[e].ijk))));
        }
      }
    }
  }

  public: void apply_forcings(int e, int n, quantity<si::time, real_t> dt)
  {
    // sanity checks
    assert(!setup->eqsys->is_homogeneous(e));
    assert(isfinite(sum((psi[e][n])(*ijk))));
    assert(isfinite(sum((rhs_R[e])(*ijk))));

    // treating explicitely the nonlinear terms
    (psi[e][n])(*ijk) += dt / si::seconds * (rhs_R[e])(*ijk);

    // treating implicitly the linear terms
    real_t C = real_t(0);
    for (int i = 0; i < setup->eqsys->var_n_rhs_terms(e); ++i)
      C += setup->eqsys->var_rhs_term(e, i).implicit_part(dt); // TODO: document that dt -> dt/2
    if (C != real_t(0))
      (psi[e][n])(*ijk) /= (real_t(1) - dt / si::seconds * C);
  }

  public: void cycle_arrays(const int e, const int n) // TODO: n unneeded?
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

  private: mtx::arr<real_t> *stash = NULL;
  private: bool stash_empty = true;
  public: void stash_cycle(int e, int n)
  {
    assert(setup->eqsys->var_dynamic(e));
    if (stash == NULL) stash = new mtx::arr<real_t>(psi[e][n].ijk); // TODO: unique_ptr
    if (stash_empty)
      *stash = psi[e][n];
    else 
      psi[e][n] = *stash;
    stash_empty = !stash_empty;
  }

  public: void integ_loop() 
  {
    assert(false); // TODO: this should not be needed with proper class hierarchy
  }
};
#endif
