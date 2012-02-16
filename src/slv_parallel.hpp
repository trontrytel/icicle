/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011 - February 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
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
  private: unique_ptr<slv_serial<real_t> > *slvs; // TODO: vector_ptr!
  private: adv<real_t> *fllbck, *advsch;
  private: stp<real_t> *setup; 
 
  public: slv_parallel(stp<real_t> *setup, out<real_t> *output,
    int i_min, int i_max,  
    int j_min, int j_max,  
    int k_min, int k_max, 
    int nsd)
    : nsd(nsd), i_min(i_min), fllbck(setup->fllbck), advsch(setup->advsch), setup(setup)
  {
    // subdomain length
    int nxl = (i_max - i_min + 1);
    nxs = nxl / nsd;
    if (nxs != ((1.*nxl) / (1.*nsd))) 
      error_macro("nxl/nsd must be an integer value (" << nxl << "/" << nsd << " given)")

    // nx=1 turns off any advection in x hence it may not be used in parallel setup
    if (nxs == 1 && nsd != 1)
      error_macro("subdomains of 1-element length not supported")

    // serial solver allocation
    slvs = new unique_ptr<slv_serial<real_t> >[nsd];
    for (int sd=0; sd < nsd; ++sd) 
      slvs[sd].reset(new slv_serial<real_t>(setup, output,
        i_min + sd * nxs, i_min + (sd + 1) * nxs - 1,
        j_min           , j_max                     , 
        k_min           , k_max                     
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
    delete[] slvs; // TODO: change to vector_ptr
  }

  private: virtual void barrier() = 0;

  private: void record(int sd, int n, unsigned long r)
  {
    barrier();
    if (sd == 0) // for shared mem the 0-th thread takes controll here
      for (int sdi=0; sdi < nsd; ++sdi) slvs[sdi]->record(n, r);
  }

  public: void integ_loop_sd(int sd) 
  {
    assert(fllbck == NULL || fllbck->num_steps() == 1);

    int n = 0;
    for (unsigned long t = 0; t < setup->nt; ++t) 
    {   
      adv<real_t> *a;
      bool fallback = this->choose_an(&a, &n, t, advsch, fllbck);

      if (t % setup->nout == 0) record(sd, n, t / setup->nout);
      for (int e = 0; e < setup->equations->n_vars(); ++e)
      {
        bool stash = // TODO: dependance on adv->time_leves... (does it work for leapfrog???)
          !setup->velocity->is_constant() && 
          a->num_steps() > 1 &&
          setup->equations->var_dynamic(e);
#  ifndef NDEBUG
        // TODO: fill caches (and perhaps also psi[n+1]) with NaNs
#  endif
        if (stash) slvs[sd]->stash_cycle(e, n); 
        for (int s = 1; s <= a->num_steps(); ++s)
        {   
          barrier(); 
          slvs[sd]->fill_halos(e, n);
          barrier(); 
          slvs[sd]->advect(e, n, s, a); 
          if (s != a->num_steps()) slvs[sd]->cycle_arrays(e, n); 
        } // s - steps
        if (stash) slvs[sd]->stash_cycle(e, n); 
      } // e - equations

      if (!setup->velocity->is_constant() && t != (setup->nt - 1)) 
      {
        barrier(); 
        // TODO: these fill haloes are repeadet in the next loop pass :(
        for (int e = 0; e < setup->equations->n_vars(); ++e)
          if (setup->equations->var_dynamic(e))
            slvs[sd]->fill_halos(e, n + 1); 
        barrier(); 
        slvs[sd]->update_courants(n + 1);
        // TODO: check if not exceeding Courant safety limits
      }

      for (int e = 0; e < setup->equations->n_vars(); ++e)
        if (!fallback) slvs[sd]->cycle_arrays(e, n);

    } // t - time
    record(sd, n, setup->nt / setup->nout);
  }

  // the two below are for MPI/fork + threads/OpenMP nested parallelisations
  public: typename mtx::arr<real_t>::type data(int e, int n, const mtx::idx &idx)
  { 
    int sd = (idx.lbound(0) - i_min) / nxs;
    assert(sd == 0 || sd == nsd - 1);
    assert((idx.ubound(0) - i_min) / nxs == sd);
    return slvs[sd]->data(e, n, idx); 
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
