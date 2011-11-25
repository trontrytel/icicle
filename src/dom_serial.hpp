/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef DOM_SERIAL_HPP
#  define DOM_SERIAL_HPP

#  include "dom.hpp"
#  include "adv.hpp"
#  include "out.hpp"
#  include "vel.hpp"
#  include "arr.hpp"
#  include "grd.hpp"

template <class unit, typename real_t>
class dom_serial : public dom<unit, real_t>
{
  private: arr<unit, real_t> *psi;
  private: Array<quantity<unit, real_t>, 3> **psi_ijk, **psi_jki, **psi_kij;
  private: adv<unit, real_t> *fllbck, *advsch;
  private: Range *i, *j, *k;
  private: out<unit, real_t> *output;
  private: int nx, ny, nz, xhalo, yhalo, zhalo;
  private: arr<si::dimensionless, real_t> *Cx, *Cy, *Cz;

  public: dom_serial(adv<unit, real_t> *fllbck, adv<unit, real_t> *advsch, 
    out<unit, real_t> *output, vel<real_t> *velocity,
    int i_min, int i_max, int nx,
    int j_min, int j_max, int ny,
    int k_min, int k_max, int nz, 
    grd<real_t> *grid,
    quantity<si::time, real_t> dt // TODO: dt should not be needed here!
  )
    : fllbck(fllbck), advsch(advsch), output(output), nx(nx), ny(ny), nz(nz)
  {
    // memory allocation
    {
      int halo = (advsch->stencil_extent() - 1) / 2;
      xhalo = halo; // in extreme cases paralellisaion may make i_max = i_min
      yhalo = (j_max != j_min ? halo : 0);
      zhalo = (k_max != k_min ? halo : 0);
    }
    psi_ijk = new Array<quantity<unit, real_t>, 3>*[advsch->time_levels()]; 
    psi_jki = new Array<quantity<unit, real_t>, 3>*[advsch->time_levels()]; 
    psi_kij = new Array<quantity<unit, real_t>, 3>*[advsch->time_levels()]; 
    for (int n=0; n < advsch->time_levels(); ++n) 
    {
      psi = new arr<unit, real_t>(
        Range(i_min - xhalo, i_max + xhalo),
        Range(j_min - yhalo, j_max + yhalo),
        Range(k_min - zhalo, k_max + zhalo)
      );
      psi_ijk[n] = &psi->ijk();
      psi_jki[n] = &psi->jki();
      psi_kij[n] = &psi->kij();
    }
    i = new Range(i_min, i_max);
    j = new Range(j_min, j_max);
    k = new Range(k_min, k_max);

    // initial condition (FIXME! - now it's always: |...... )
    //for (int i_int = i->first(); i_int <= i->last(); ++i_int)
    //  for (int j_int = j->first(); j_int <= j->last(); ++j_int)
    //    for (int k_int = k->first(); k_int <= k->last(); ++k_int)
    //      if (i_int == k_int) (*psi_ijk[0])(i_int, j_int, k_int) = 1;
    if (i_min == 0) (*psi_ijk[0])(Range(0),0,0) = 1;

    // periodic boundary
    this->hook_neighbours(this, this);
 
    // velocity fields
    {
      Range 
        ii = grid->rng_vctr(i->first(), i->last()),
        jj = (j->first() != j->last())
          ? grid->rng_vctr(j->first(), j->last())
          : Range(j->first(), j->last()),
        kk = (k->first() != k->last())
          ? grid->rng_vctr(k->first(), k->last())
          : Range(k->first(), k->last());

      Cx = new arr<si::dimensionless, real_t>(ii, jj, kk);
      Cy = new arr<si::dimensionless, real_t>(ii, jj, kk);
      Cz = new arr<si::dimensionless, real_t>(ii, jj, kk);

      grid->populate_courant_fields(ii, jj, kk, &Cx->ijk(), &Cy->ijk(), &Cz->ijk(), velocity, dt);
    }
  }

  public: ~dom_serial()
  {
    for (int n=0; n < advsch->time_levels(); ++n) delete psi_jki[n];
    delete[] psi_jki;
    for (int n=0; n < advsch->time_levels(); ++n) delete psi_kij[n];
    delete[] psi_kij;
    for (int n=0; n < advsch->time_levels(); ++n) delete psi_ijk[n];
    delete[] psi_ijk;
    delete psi;
    delete i;
    delete Cx; delete Cy; delete Cz;
  }

  public: void record(const int n, const unsigned long t)
  {
    output->record(psi_ijk, n, *i, *j, *k, t);
  }

  private: quantity<unit, real_t> data(int n, int i, int j, int k) 
  { 
    return (*psi_ijk[n])(i, j, k);
  }

  public: void fill_halos(int n)
  {
    // left halo
    for (int i_int = i->first() - xhalo; i_int < i->first(); ++i_int)
      for (int j_int = j->first(); j_int <= j->last(); ++j_int)
        for (int k_int = k->first(); k_int <= k->last(); ++k_int)
          (*psi_ijk[n])(i_int, j_int, k_int) = 
            this->left_nghbr_data(n, (i_int + nx) % nx, j_int, k_int);

    // rght halo
    for (int i_int = i->last() + 1; i_int <= i->last() + xhalo; ++i_int)
      for (int j_int = j->first(); j_int <= j->last(); ++j_int)
        for (int k_int = k->first(); k_int <= k->last(); ++k_int)
          (*psi_ijk[n])(i_int, j_int, k_int) = 
            this->rght_nghbr_data(n, (i_int + nx) % nx, j_int, k_int);
  }

  public: void advect(adv<unit, real_t> *a, int n, int s, quantity<si::time, real_t> dt)
  {
    // op_1D() uses the -= operator so the first assignment happens here
    *psi_ijk[n+1] = *psi_ijk[0]; 

    if (true) // in extreme cases paralellisaion may make i->first() = i->last()
      a->op_1D(psi_ijk, *i, *j, *k, n, s, Cx->ijk(), Cy->ijk(), Cz->ijk()); // X
    if (j->first() != j->last()) 
      a->op_1D(psi_jki, *j, *k, *i, n, s, Cy->jki(), Cz->jki(), Cx->jki()); // Y
    if (k->first() != k->last()) 
      a->op_1D(psi_kij, *k, *i, *j, n, s, Cz->kij(), Cx->kij(), Cy->kij()); // Z
  }

  public: void cycle_arrays(const int n)
  {
    switch (advsch->time_levels())
    {
      case 2: // e.g. upstream, mpdata
        cycleArrays(*psi_ijk[n], *psi_ijk[n+1]); 
        cycleArrays(*psi_jki[n], *psi_jki[n+1]); 
        cycleArrays(*psi_kij[n], *psi_kij[n+1]); 
        break;
      case 3: // e.g. leapfrog
        cycleArrays(*psi_ijk[n-1], *psi_ijk[n], *psi_ijk[n+1]); 
        cycleArrays(*psi_jki[n-1], *psi_jki[n], *psi_jki[n+1]); 
        cycleArrays(*psi_kij[n-1], *psi_kij[n], *psi_kij[n+1]); 
        break;
      default: assert(false);
    }
  }

  public: void integ_loop(unsigned long nt, quantity<si::time, real_t> dt) 
  {
    int n;
    for (unsigned long t = 0; t < nt; ++t)
    {   
      adv<unit, real_t> *a;
      bool fallback = this->choose_an(&a, &n, t, advsch, fllbck);    
      record(n, t);
      for (int s = 1; s <= a->num_steps(); ++s)
      {   
        fill_halos(n);
        advect(a, n, s, dt); 
        if (!fallback) cycle_arrays(n);
      } // s
    } // t
    record(n, nt);
  }
};
#endif
