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
  private: ptr_vector<slv_serial<real_t>> slvs; 
  private: const stp<real_t> &setup; 
 
  public: slv_parallel(
    const stp<real_t> &setup, 
    out<real_t> &output,
    int i_min, int i_max,  
    int j_min, int j_max,  
    int k_min, int k_max, 
    int nsd)
    : nsd(nsd), i_min(i_min), setup(setup)
  {
    // subdomain length
    int nxl = (i_max - i_min + 1);
    nxs = nxl / nsd;
    if (nxs != ((1.*nxl) / (1.*nsd))) 
      error_macro("nxl/nsd must be an integer value (" << nxl << "/" << nsd << " given)")

    // nx=1 turns off any advection in x hence it may not be used in parallel setup
    if (nxs == 1 && nsd != 1)
      error_macro("subdomains of 1-element length not supported")

    // serial solver allocation (includes loading of initial condition)
    for (int sd=0; sd < nsd; ++sd) 
      slvs.push_back(new slv_serial<real_t>(setup, output,
        i_min + sd * nxs, i_min + (sd + 1) * nxs - 1,
        j_min           , j_max                     , 
        k_min           , k_max                     
      ));

    // periodic boundary over prallel domain
    for (int sd=0; sd < nsd; ++sd) 
    {   
      int l = (nsd + sd - 1) % nsd; 
      int r = (nsd + sd + 1) % nsd;  
      slvs[sd].hook_neighbour(slv<real_t>::left, &slvs[l]);
      slvs[sd].hook_neighbour(slv<real_t>::rght, &slvs[r]);
    }
  }

  private: virtual void barrier() = 0;

  private: void record(int sd, int n, unsigned long r)
  {
    barrier();
    if (sd == 0) // for shared mem the 0-th thread takes controll here
      for (int sdi=0; sdi < nsd; ++sdi) slvs[sdi].record(n, r);
  }

  private: void fill_halos(int sd, int e, int n)
  {
    barrier();
    slvs[sd].fill_halos(e, n); 
    barrier(); // TODO: not needed with distmem?
  }

  private: void advect(int sd, int e, int n)
  {
    for (int s = 1; s <= setup.advsch->num_steps(); ++s)
    {   
      if (s != 1) fill_halos(sd, e, n);
      slvs[sd].advect(e, n, s); 
      if (s != setup.advsch->num_steps()) slvs[sd].cycle_arrays(e, n); 
    }
  }

  public: void integ_loop_sd(int sd) 
  {
    vector<bool> homo(setup.eqsys->n_vars());
    bool allhomo = true;
    for (int e = 0; e < setup.eqsys->n_vars(); ++e)
    {
      homo[e] = setup.eqsys->is_homogeneous(e);
      if (!homo[e]) allhomo = false;
    }
    int n = 0;

    // adjustments for the initial condition
    slvs[sd].apply_adjustments(n, setup.dt);

    for (int e = 0; e < setup.eqsys->n_vars(); ++e) 
      fill_halos(sd, e, n); 

    // first guess for velocity fields assuming constant momenta
    if (!setup.velocity->is_constant())
      slvs[sd].copy(n, n + 1); // TODO: only "dynamic" variables

    // first guess for forcing terms at t=0
    if (!allhomo) slvs[sd].update_forcings(n, 0 * si::seconds);

    // time stepping
    for (unsigned long t = 0; t < setup.nt; ++t) 
    {   
      if (sd == 0) cerr << "t/dt=" << t << " (t=" << real_t(t) * setup.dt << ")" << endl;

      assert(setup.advsch->time_levels() <= 3); // TODO: support for other values
      bool fallback = (t == 0 && setup.advsch->time_levels() == 3);
      n = fallback ? 0 : setup.advsch->time_levels() - 2;

      if (t % setup.nout == 0) record(sd, n, t / setup.nout);

      int last_group = -1;
      for (int e = 0; e < setup.eqsys->n_vars(); ++e)
      {
        // filling halos
        for (int e = 0; e < setup.eqsys->n_vars(); ++e) 
          fill_halos(sd, e, n); 

        // updating velocity field
        {
          int group = setup.eqsys->group(e);
          if (!setup.velocity->is_constant() && group != last_group)
          {
            slvs[sd].update_courants(group, n + 1, n); // TODO: vector of n-s
            last_group = group;
          }
        }

        // TODO: dependance on adv->time_leves... (does it work for leapfrog???)
        bool stash = (setup.advsch->num_steps() > 1 && setup.eqsys->var_dynamic(e))
          || (!setup.velocity->is_constant() && !homo[e]);

        slvs[sd].fill_with_nans(e, n + 1); // TODO: fill caches with NaNs?
        if (stash) slvs[sd].stash_cycle(e, n); 
        if (!homo[e]) slvs[sd].apply_forcings(e, n, real_t(.5) * setup.dt);
        advect(sd, e, n);
        if (stash) slvs[sd].stash_cycle(e, n); 
      } // e - equations

      // apply post-advection, pre-rhs adjustments
      slvs[sd].apply_adjustments(n + 1, setup.dt);

      // assuming that computation of forcings relies on the outcome of advection of ALL variables
      if (!allhomo)
      {
        for (int e = 0; e < setup.eqsys->n_vars(); ++e) 
          fill_halos(sd, e, n + 1); 
        slvs[sd].update_forcings(n + 1, real_t(t + .5) * setup.dt);
        for (int e = 0; e < setup.eqsys->n_vars(); ++e)
          if (!homo[e]) slvs[sd].apply_forcings(e, n + 1, real_t(.5) * setup.dt);
      } 

      if (!fallback) 
        for (int e = 0; e < setup.eqsys->n_vars(); ++e)
          slvs[sd].cycle_arrays(e, n);
    } // t - time
    record(sd, n, setup.nt / setup.nout);
  }

  // the two below are for MPI/fork + threads/OpenMP nested parallelisations
  public: typename mtx::arr<real_t>::type data(int e, int n, const mtx::idx &idx)
  { 
    int sd = (idx.lbound(0) - i_min) / nxs;
    assert(sd == 0 || sd == nsd - 1);
    assert((idx.ubound(0) - i_min) / nxs == sd);
    return slvs[sd].data(e, n, idx); 
  }

  public: void hook_neighbour(int s, slv<real_t> *n) 
  {
    switch (s)
    {
      case slv<real_t>::left: slvs[0    ].hook_neighbour(s, n); break;
      case slv<real_t>::rght: slvs[nsd-1].hook_neighbour(s, n); break;
      default: assert(false);
    }
  }
};
#endif
