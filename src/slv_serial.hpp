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

template <typename real_t>
class slv_serial : public slv<real_t>
{
  private: adv<real_t> *fllbck, *advsch;
  private: auto_ptr<Range> i, j, k;
  private: out<real_t> *output;
  private: int nx, ny, nz, halo;//, halo, zhalo;

  private: auto_ptr<Array<real_t, 3> > Cx, Cy, Cz;

  private: auto_ptr<Array<real_t, 3> > *psi_guard;
  private: Array<real_t, 3> **psi;

  private: auto_ptr<Array<real_t, 3> > *tmp_s_guard;
  private: Array<real_t, 3> **tmp_s;
  private: auto_ptr<Array<real_t, 3> > *tmp_v_guard;
  private: Array<real_t, 3> **tmp_v;

  public: slv_serial(stp<real_t> *setup,
    int i_min, int i_max, int nx,
    int j_min, int j_max, int ny,
    int k_min, int k_max, int nz, 
    quantity<si::time, real_t> dt // TODO: dt should not be needed here!
  )
    : fllbck(setup->fllbck), advsch(setup->advsch), output(setup->output), nx(nx), ny(ny), nz(nz)
  {
    halo = (advsch->stencil_extent() - 1) / 2;

    // sanity checks
    {
      int len;
      if (halo > (len = i_max - i_min + 1)) 
        error_macro("halo length (" << halo << ") may not exceed domain extent in X (" << len << ")")
      if (halo > (len = j_max - j_min + 1))
        error_macro("halo length (" << halo << ") may not exceed domain extent in X (" << len << ")")
      if (halo > (len = k_max - k_min + 1))
        error_macro("halo length (" << halo << ") may not exceed domain extent in X (" << len << ")")
    }

    // memory allocation
    psi_guard = new auto_ptr<Array<real_t, 3> >[advsch->time_levels()];
    psi = new Array<real_t, 3>*[advsch->time_levels()];
    for (int n=0; n < advsch->time_levels(); ++n) 
    {
      psi_guard[n].reset(new Array<real_t, 3>(
        setup->grid->rng_sclr(i_min, i_max, halo),
        setup->grid->rng_sclr(j_min, j_max, halo),
        setup->grid->rng_sclr(k_min, k_max, halo)
      ));
      psi[n] = psi_guard[n].get();
    }

    // caches
    assert(fllbck == NULL || fllbck->num_sclr_caches() == 0);
    assert(fllbck == NULL || fllbck->num_vctr_caches() == 0);

    {
      int cnt = advsch->num_vctr_caches();
      tmp_v_guard = new auto_ptr<Array<real_t, 3> >[cnt];
      tmp_v = new Array<real_t, 3>*[cnt];
      for (int n=0; n < cnt; ++n)
      {
        tmp_v_guard[n].reset(new Array<real_t, 3>(
          setup->grid->rng_vctr(i_min, i_max),
          setup->grid->rng_vctr(j_min, j_max),
          setup->grid->rng_vctr(k_min, k_max)
        ));
        tmp_v[n] = tmp_v_guard[n].get();
      }
    }
    {
      int cnt = advsch->num_sclr_caches();
      tmp_s_guard = new auto_ptr<Array<real_t, 3> >[cnt];
      tmp_s = new Array<real_t, 3>*[cnt];
      for (int n=0; n < cnt; ++n)
      {
        tmp_s_guard[n].reset(new Array<real_t, 3>(
          setup->grid->rng_vctr(i_min, i_max),
          setup->grid->rng_vctr(j_min, j_max),
          setup->grid->rng_vctr(k_min, k_max)
        ));
        tmp_s[n] = tmp_s_guard[n].get();
      }
    }

    // indices
    i.reset(new Range(i_min, i_max));
    j.reset(new Range(j_min, j_max));
    k.reset(new Range(k_min, k_max));

    // initial condition
    setup->grid->populate_scalar_field(*i, *j, *k, psi[0], setup->intcond);

    // periodic boundary in all directions
    for (int s=this->first; s <= this->last; ++s) 
      this->hook_neighbour(s, this);
 
    // velocity fields
    {
      Range 
        ii = setup->grid->rng_vctr(i->first(), i->last()),
        jj = setup->grid->rng_vctr(j->first(), j->last()),
        kk = setup->grid->rng_vctr(k->first(), k->last());
      //TODO dla nx=3 ny=3 Cx powinno być 4x3 a nie 4x4
      Cx.reset(new Array<real_t, 3>(ii, jj, kk));
      Cy.reset(new Array<real_t, 3>(ii, jj, kk));
      Cz.reset(new Array<real_t, 3>(ii, jj, kk));

      setup->grid->populate_courant_fields(ii, jj, kk, Cx.get(), Cy.get(), Cz.get(), setup->velocity, dt);
    }
  }

  public: ~slv_serial()
  {
    delete[] psi;
    delete[] psi_guard;
    //delete[] tmp_{s,v}; TODO!!!
    //delete[] tmp_{s,v}_guard; TODO!!!
  }

  public: void record(const int n, const unsigned long t)
  {
    output->record(psi[n], *i, *j, *k, t);
  }

  public: Array<real_t, 3> data(int n, const RectDomain<3> &idx)
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
  void fill_halos_helper(Array<real_t, 3> *psi[], int nghbr, int n, 
    int i_min, int i_max, const Range &j, const Range &k, int mod
  )
  {
    (*psi[n])(idx(Range(i_min, i_max), j, k)) =
      this->nghbr_data(nghbr, n, idx(Range((i_min + mod) % mod, (i_max + mod) % mod), j, k));
  }

  public: void advect(adv<real_t> *a, int n, int s, quantity<si::time, real_t> dt)
  {
    a->op3D(psi, tmp_s, tmp_v, *i, *j, *k, n, s, *Cx, *Cy, *Cz);
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
