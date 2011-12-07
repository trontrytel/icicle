/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef SLV_PARALLEL_HPP
#  define SLV_PARALLEL_HPP

#  include "slv_serial.hpp"
#  include "out.hpp"

template <typename real_t>
class slv_parallel : public slv<real_t>
{
  private: int nsd; // number of subdomains within this solver
  private: int nxs; // number of points per subdomain
  private: int i_min; // i_min for this set of subdomains
  private: auto_ptr<slv_serial<real_t> > *slvs;
  private: adv<real_t> *fllbck, *advsch;
  
  public: slv_parallel(stp<real_t> *setup,
    int i_min, int i_max, int nx, 
    int j_min, int j_max, int ny, 
    int k_min, int k_max, int nz,
    quantity<si::time, real_t> dt,
    int nsd)
    : nsd(nsd), i_min(i_min), fllbck(setup->fllbck), advsch(setup->advsch)
  {
    // subdomain length
    int nxl = (i_max - i_min + 1);
    nxs = nxl / nsd;
    if (nxs != ((1.*nxl) / (1.*nsd))) 
      error_macro("nxl/nsd must be an integer value (" << nxl << "/" << nsd << " given)")

    // serial solver allocation (TODO: there could be just one psi for all subdomains...)
    slvs = new auto_ptr<slv_serial<real_t> >[nsd];
    for (int sd=0; sd < nsd; ++sd) 
      slvs[sd].reset(new slv_serial<real_t>(setup,
        i_min + sd * nxs, i_min + (sd + 1) * nxs - 1, nx,
        j_min           , j_max                     , ny, 
        k_min           , k_max                     , nz,
        dt
      ));

    // periodic boundary over prallel domain
    for (int sd=0; sd < nsd; ++sd) 
    {   
      int l = (nsd + sd - 1) % nsd; 
      int r = (nsd + sd + 1) % nsd;  
      slvs[sd]->hook_neighbour(slv<real_t>::left, slvs[l].get());
      slvs[sd]->hook_neighbour(slv<real_t>::rght, slvs[r].get());
    }
  }

  public: ~slv_parallel()
  {
    delete[] slvs; // vector is not used as it might not work with auto_ptr
  }

  private: virtual void barrier() = 0;

  public: void integ_loop_sd(unsigned long nt, quantity<si::time, real_t> dt, int sd) 
  {
    int n = 0;
    for (unsigned long t = 0; t < nt; ++t)
    {   
      adv<real_t> *a;
      bool fallback = this->choose_an(&a, &n, t, advsch, fllbck);
      barrier();
      if (sd == 0) for (int sdi=0; sdi < nsd; ++sdi) slvs[sdi]->record(n, t);
      for (int s = 1; s <= a->num_steps(); ++s)
      {   
        barrier();
        slvs[sd]->fill_halos(n);
        barrier();
        slvs[sd]->advect(a, n, s, dt); 
        if (!fallback) slvs[sd]->cycle_arrays(n);
      } // s
    } // t
    barrier();
    if (sd == 0) for (int sd=0; sd < nsd; ++sd) slvs[sd]->record(n, nt);
  }

  // the two below are for MPI/fork + threads/OpenMP nested parallelisations
  public: Array<real_t, 3> data(int n, 
    const Range &i, const Range &j, const Range &k) 
  { 
    int sd = (i.first() - i_min) / nxs;
    assert(sd == 0 || sd == nsd - 1);
    assert((i.last() - i_min) / nxs == sd);
    return slvs[sd]->data(n, i, j, k);
  }

  public: void hook_neighbour(int s, slv<real_t> *n) 
  {
    switch (s)
    {
      case slv<real_t>::left: slvs[0    ]->hook_neighbour(s, n); break;
      case slv<real_t>::rght: slvs[nsd-1]->hook_neighbour(s, n); break;
      default: assert(false);
    }
  }
};
#endif