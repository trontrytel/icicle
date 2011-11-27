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

template <class unit, typename real_t>
class slv_parallel : public slv<unit, real_t>
{
  private: int nsd;
  private: auto_ptr<slv_serial<unit, real_t> > *slvs;
  private: adv<unit, real_t> *fllbck, *advsch;
  
  public: slv_parallel(adv<unit, real_t> *fllbck, adv<unit, real_t> *advsch, 
    out<unit, real_t> *output, vel<real_t> *velocity, 
    int nx, int ny, int nz,
    grd<real_t> *grid,
    quantity<si::time, real_t> dt,
    int nsd)
    : nsd(nsd), fllbck(fllbck), advsch(advsch)
  {
    // subdomain length
    int nxs = nx / nsd;
    if (nxs != ((1.*nx) / (1.*nsd))) 
      error_macro("nx/nk must be an integer value (" << nx << "/" << nsd << " given)")

    // serial solver allocation
    slvs = new auto_ptr<slv_serial<unit, real_t> >[nsd];
    for (int sd=0; sd < nsd; ++sd) 
      slvs[sd].reset(new slv_serial<unit, real_t>(fllbck, advsch, output, velocity,
        sd * nxs, (sd + 1) * nxs - 1, nx,
        0,        ny - 1,             ny, 
        0,        nz - 1,             nz, 
        grid,
        dt
      ));

    // periodic boundary over prallel domain
    for (int sd=0; sd < nsd; ++sd) 
    {   
      int l = (nsd + sd - 1) % nsd; 
      int r = (nsd + sd + 1) % nsd;  
      slvs[sd]->hook_neighbours(slvs[l].get(), slvs[r].get());
    }
  }

  public: ~slv_parallel()
  {
    delete[] slvs;
  }

  private: virtual void barrier() = 0;

  public: void integ_loop_sd(unsigned long nt, quantity<si::time, real_t> dt, int sd) 
  {
    int n = 0;
    for (unsigned long t = 0; t < nt; ++t)
    {   
      adv<unit, real_t> *a;
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

  // e.g. for MPI + OpenMP nesting (FIXME)
  private: quantity<unit, real_t> data(int n, int i, int j, int k) { throw; }
};
#endif
