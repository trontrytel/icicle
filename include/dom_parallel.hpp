/** @file
 *  @author Sylwester Arabas (University of Warsaw) <slayoo@igf.fuw.edu.pl>
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef DOM_PARALLEL_HPP
#  define DOM_PARALLEL_HPP

#  include "dom.hpp"
#  include "adv.hpp"
#  include "out.hpp"

template <class unit, typename real_t>
class dom_parallel : public dom<unit, real_t>
{
  private: int nsd;
  private: dom_serial<unit, real_t> **doms;
  private: adv<unit, real_t> *fllbck, *advsch;
  
  public: dom_parallel(adv<unit, real_t> *fllbck, adv<unit, real_t> *advsch, 
    out<unit, real_t> *output, int nx, int ny, int nz, int nsd)
    : nsd(nsd), fllbck(fllbck), advsch(advsch)
  {
    // subdomain length
    int nxs = nx / nsd;
    if (nxs != ((1.*nx) / (1.*nsd))) 
      error_macro("nx/nk must be an integer value (" << nx << "/" << nsd << " given)")

    // domain allocation
    doms = new dom_serial<unit, real_t>*[nsd];
    for (int sd=0; sd < nsd; ++sd) 
      doms[sd] = new dom_serial<unit, real_t>(fllbck, advsch, output, 
        sd * nxs, (sd + 1) * nxs - 1, nx,
        0,        ny - 1,             ny,
        0,        nz - 1,             nz
      );

    // periodic boundary over prallel domain
    for (int sd=0; sd < nsd; ++sd) 
    {   
      int l = (nsd + sd - 1) % nsd; 
      int r = (nsd + sd + 1) % nsd;  
      doms[sd]->hook_neighbours(doms[l], doms[r]);
    }
  }

  public: ~dom_parallel()
  {
    for (int sd=0; sd < nsd; ++sd) delete doms[nsd];
    delete[] doms;
  }

  private: virtual void barrier() = 0;

  public: void integ_loop_sd(unsigned long nt, 
    quantity<si::dimensionless, real_t> &Cx, 
    quantity<si::dimensionless, real_t> &Cy, 
    quantity<si::dimensionless, real_t> &Cz, 
    int sd) 
  {
    int n = 0;
    for (unsigned long t = 0; t < nt; ++t)
    {   
      adv<unit, real_t> *a;
      bool fallback = choose_an(&a, &n, t, advsch, fllbck);
      barrier();
      if (sd == 0) for (int sdi=0; sdi < nsd; ++sdi) doms[sdi]->record(n, t);
      for (int s = 1; s <= a->num_steps(); ++s)
      {   
        barrier();
        doms[sd]->fill_halos(n);
        barrier();
        doms[sd]->advect(a, n, Cx, Cy, Cz, s); 
        if (!fallback) doms[sd]->cycle_arrays(n);
      } // s
    } // t
    barrier();
    if (sd == 0) for (int sd=0; sd < nsd; ++sd) doms[sd]->record(n, nt);
  }

  // e.g. for MPI + OpenMP nesting (FIXME)
  private: quantity<unit, real_t> data(int n, int i, int j, int k) { throw; }
};
#endif
