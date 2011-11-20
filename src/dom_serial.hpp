/** @file
 *  @author Sylwester Arabas (University of Warsaw) <slayoo@igf.fuw.edu.pl>
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

template <class unit, typename real_t>
class dom_serial : public dom<unit, real_t>
{
  private: Array<quantity<unit, real_t>, 3> 
    **psi_ijk, **psi_jki, **psi_kij;
  private: adv<unit, real_t> *fllbck, *advsch;
  private: Range *i, *j, *k;
  private: out<unit, real_t> *output;
  private: int nx, ny, nz, xhalo, yhalo, zhalo;
  private: Array<quantity<si::dimensionless, real_t>, 3> 
    *Cx_ijk, *Cy_ijk, *Cz_ijk,
    *Cx_jki, *Cy_jki, *Cz_jki,
    *Cx_kij, *Cy_kij, *Cz_kij;
  private: quantity<si::length, real_t> dx, dy, dz;

  // TODO: move to an 'arr' class that would have ijk(), jki() and kij() methods
  private: Array<quantity<unit, real_t>, 3>*array_view(
    Array<quantity<unit, real_t>, 3>* a, int d1, int d2, int d3
  )
  {
    GeneralArrayStorage<3> storage;
    storage.ordering() = d1, d2, d3;
    storage.base() = a->lbound(d1), a->lbound(d2), a->lbound(d3);
    return new Array<quantity<unit, real_t>, 3>(
      a->dataFirst(), 
      TinyVector<int, 3>(a->extent(d1), a->extent(d2), a->extent(d3)), 
      storage
    );
  }

  public: dom_serial(adv<unit, real_t> *fllbck, adv<unit, real_t> *advsch, 
    out<unit, real_t> *output, vel<real_t> *velocity,
    int i_min, int i_max, int nx, quantity<si::length, real_t> dx,
    int j_min, int j_max, int ny, quantity<si::length, real_t> dy,
    int k_min, int k_max, int nz, quantity<si::length, real_t> dz,
    quantity<si::time, real_t> dt // TODO: dt should not be needed here!
  )
    : fllbck(fllbck), advsch(advsch), output(output), nx(nx), ny(ny), nz(nz), dx(dx), dy(dy), dz(dz)
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
      psi_ijk[n] = new Array<quantity<unit, real_t>, 3>(
        Range(i_min - xhalo, i_max + xhalo),
        Range(j_min - yhalo, j_max + yhalo),
        Range(k_min - zhalo, k_max + zhalo)
      ); 

      // two "views" of the array with permuted dimensions (but the data remains the same!)
      psi_jki[n] = array_view(psi_ijk[n], secondRank, thirdRank, firstRank);
      psi_kij[n] = array_view(psi_ijk[n], thirdRank, firstRank, secondRank);
    }
    i = new Range(i_min, i_max);
    j = new Range(j_min, j_max);
    k = new Range(k_min, k_max);

    // initial condition (FIXME! - now it's always: |...... )
    if (i_min == 0) (*psi_ijk[0])(Range(0),0,0) = 1;

    // periodic boundary
    this->hook_neighbours(this, this);
 
    // velocity fields
    {
      Range 
        ii = Range(i->first() - grd::m_half, i->last() + grd::p_half),
        jj = (j->first() != j->last())
          ? Range(j->first() - grd::m_half, j->last() + grd::p_half)
          : Range(j->first(), j->last()),
        kk = (k->first() != k->last())
          ? Range(k->first() - grd::m_half, k->last() + grd::p_half)
          : Range(k->first(), k->last());

      Cx_ijk = new Array<quantity<si::dimensionless, real_t>, 3>(ii, jj, kk);
      Cx_jki = array_view(Cx_ijk, secondRank, thirdRank, firstRank);
      Cx_kij = array_view(Cx_ijk, thirdRank, firstRank, secondRank);

      Cy_ijk = new Array<quantity<si::dimensionless, real_t>, 3>(ii, jj, kk);
      Cy_jki = array_view(Cy_ijk, secondRank, thirdRank, firstRank);
      Cy_kij = array_view(Cy_ijk, thirdRank, firstRank, secondRank);

      Cz_ijk = new Array<quantity<si::dimensionless, real_t>, 3>(ii, jj, kk);
      Cz_jki = array_view(Cz_ijk, secondRank, thirdRank, firstRank);
      Cz_kij = array_view(Cz_ijk, thirdRank, firstRank, secondRank);

      for (int i_int = ii.first(); i_int <= ii.last(); ++i_int)
      {
        for (int j_int = jj.first(); j_int <= jj.last(); ++j_int)
        {
          for (int k_int = kk.first(); k_int <= kk.last(); ++k_int)
          {
            quantity<si::length, real_t> 
              x = real_t(i_int - .5) * dx,
              y = real_t(j_int - .5) * dy,
              z = real_t(k_int - .5) * dz;
            (*Cx_ijk)(i_int, j_int, k_int) = velocity->u(x, y, z) * dt / dx;
            if (j->first() != j->last()) 
              (*Cy_ijk)(i_int, j_int, k_int) = velocity->v(x, y, z) * dt / dy;
            if (k->first() != k->last()) 
              (*Cz_ijk)(i_int, j_int, k_int) = velocity->w(x, y, z) * dt / dz;
          }
        }
      }
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
    delete i;
    delete Cx_ijk; delete Cy_ijk; delete Cz_ijk;
    delete Cx_jki; delete Cy_jki; delete Cz_jki;
    delete Cx_kij; delete Cy_kij; delete Cz_kij;
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
    *psi_ijk[n+1] = *psi_ijk[0]; // op_1D() uses the -= operator so the first assignment happens here
    a->op_1D(psi_ijk, *i, n, s, *Cx_ijk, *Cy_ijk, *Cz_ijk); // in extreme cases paralellisaion may make i->first() = i->last()
    if (j->first() != j->last()) a->op_1D(psi_jki, *j, n, s, *Cy_jki, *Cz_jki, *Cx_jki);
    if (k->first() != k->last()) a->op_1D(psi_kij, *k, n, s, *Cz_kij, *Cx_kij, *Cy_kij);
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
