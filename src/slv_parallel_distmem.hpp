/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef SLV_PARALLEL_DISTMEM_HPP
#  define SLV_PARALLEL_DISTMEM_HPP

#  include "slv_parallel.hpp" 

template <typename real_t, class shrdmem_class>
class slv_parallel_distmem : public shrdmem_class
{
  // <nested class>
  private: class slv_halo : public slv<real_t>
  {
    private: slv_parallel_distmem<real_t, shrdmem_class> *nghbr;
    private: int peer, cnt;
    private: auto_ptr<Array<real_t, 3> > ibuf, obuf;
    private: Range ixr, oxr, yr, zr;

    public: slv_halo(slv_parallel_distmem<real_t, shrdmem_class> *nghbr, int peer,
      const Range &ixr, const Range &oxr, const Range &yr, const Range &zr
    )
      : nghbr(nghbr), peer(peer), ixr(ixr), oxr(oxr), yr(yr), zr(zr)
    { 
      ibuf.reset(new Array<real_t, 3>(ixr, yr, zr));
      obuf.reset(new Array<real_t, 3>(oxr, yr, zr));
      cnt = ibuf->cols() * ibuf->rows() * ibuf->depth();
    }

    public: Array<real_t, 3> data(int n, const Range &i, const Range &j, const Range &k) 
    { 
      return (*ibuf)(i, j, k); 
    }
 
    public: void sync(int n)
    {
      *obuf = nghbr->data(n, oxr, yr, zr);
      nghbr->sndrcv(peer, cnt, ibuf->data(), obuf->data());
    }

    public: void integ_loop(long unsigned, quantity<si::time, real_t>) 
    { 
      assert(false); 
    }
  };
  // </nested class>

  private: int size, rank;
  private: auto_ptr<slv_halo> lhalo, rhalo;
  private: auto_ptr<Array<real_t, 3> > libuf, ribuf, lobuf, robuf;

  private: int i_min(int nx, int rank, int size) { return (rank + 0) * nx / size; }
  private: int i_max(int nx, int rank, int size) { return i_min(nx, rank + 1, size) - 1; }

  public: slv_parallel_distmem(stp<real_t> *setup,
    int nx, int ny, int nz, quantity<si::time, real_t> dt, int size, int rank)
    : shrdmem_class(setup,
      i_min(nx, rank, size), i_max(nx, rank, size), nx,
      0, ny - 1, ny,
      0, nz - 1, nz,
      dt,
      1 // FIXME: MPI+OpenMP and MPI+threads 
    ), size(size), rank(rank)
  {
    // we distinguish left and right halos just by knowing where they come from so...
    if (size <= 2) error_macro("at least three subdomains are needed for distmem parallelisations")

    // subdomain length
    int nxs = nx / size;
    if (nxs != ((1.*nx) / (1.*size))) error_macro("nx/nk must be an integer value (" << nx << "/" << size << " given)")

    // halo containers for halo domain (TODO: use some clever Blitz constructor to save memory)
    {
      int halo = (setup->advsch->stencil_extent() - 1) / 2;
      int peer_left = (size + rank - 1) % size;
      int peer_rght = (size + rank + 1) % size;

      Range ixr, oxr, yr(0, ny - 1), zr(nz - 1);

      ixr = Range( // input xr (halo)
        (i_min(nx, rank, size) - halo + nx) % nx, 
        (i_min(nx, rank, size) - 1    + nx) % nx  
      );
      oxr = Range( // output xr (edge)
        i_min(nx, rank, size), 
        i_min(nx, rank, size) + halo - 1
      );
      lhalo.reset(new slv_halo(this, peer_left, ixr, oxr, yr, zr));

      ixr = Range( // input xr (halo)
        (i_max(nx, rank, size) + 1    + nx) % nx, 
        (i_max(nx, rank, size) + halo + nx) % nx
      );
      oxr = Range( // output xr (edge)
        i_max(nx, rank, size) - halo + 1,
        i_max(nx, rank, size)  
      );
      rhalo.reset(new slv_halo(this, peer_rght, ixr, oxr, yr, zr));
    }
    this->hook_neighbour(slv<real_t>::left, lhalo.get());
    this->hook_neighbour(slv<real_t>::rght, rhalo.get());
  }

  private: virtual void distmem_barrier() = 0;

  private: void barrier() 
  {
    shrdmem_class::barrier();
    distmem_barrier();
  }

  private: virtual void sndrcv(int peer, int cnt, real_t *ibuf, real_t *obuf) = 0;
};

#endif
