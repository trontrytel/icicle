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
    private: auto_ptr<arr<real_t> > ibuf, obuf;
    private: rng ixr, oxr, yr, zr;

    public: slv_halo(slv_parallel_distmem<real_t, shrdmem_class> *nghbr, int peer,
      const rng &ixr, const rng &oxr, const rng &yr, const rng &zr
    )
      : nghbr(nghbr), peer(peer), ixr(ixr), oxr(oxr), yr(yr), zr(zr)
    { 
      ibuf.reset(new arr<real_t>(ixr, yr, zr));
      obuf.reset(new arr<real_t>(oxr, yr, zr));
      cnt = ibuf->cols() * ibuf->rows() * ibuf->depth();
    }

    public: typename arr<real_t>::arr_ret data(int n, const idx &idx)
    { 
      return (*ibuf)(idx); 
    }
 
    public: void sync(int n)
    {
      (*obuf)(ixr, yr, zr) = nghbr->data(n, idx_ijk(oxr, yr, zr));
      nghbr->sndrcv(peer, cnt, ibuf->data(), obuf->data());
    }

    public: void integ_loop() 
    { 
      assert(false); 
    }
  };
  // </nested class>

  private: int size, rank;
  private: auto_ptr<slv_halo> lhalo, rhalo;
  private: auto_ptr<arr<real_t> > libuf, ribuf, lobuf, robuf;

  private: int i_min(int nx, int rank, int size) { return (rank + 0) * nx / size; }
  private: int i_max(int nx, int rank, int size) { return i_min(nx, rank + 1, size) - 1; }

  public: slv_parallel_distmem(stp<real_t> *setup, out<real_t> *output,
    int size, int rank)
    : shrdmem_class(setup, output,
      i_min(setup->nx, rank, size), i_max(setup->nx, rank, size), 
      0, setup->ny - 1,
      0, setup->nz - 1,
      1 // FIXME: MPI+OpenMP and MPI+threads 
    ), size(size), rank(rank)
  {
    // we distinguish left and right halos just by knowing where they come from so...
    if (size <= 2) error_macro("at least three subdomains are needed for distmem parallelisations")

    // subdomain length
    int nxs = setup->nx / size;
    if (nxs != ((1.*setup->nx) / (1.*size))) error_macro("nx/nk must be an integer value (" << setup->nx << "/" << size << " given)")

    // halo containers for halo domain (TODO: use some clever Blitz constructor to save memory)
    {
      int halo = (setup->advsch->stencil_extent() - 1) / 2;
      int peer_left = (size + rank - 1) % size;
      int peer_rght = (size + rank + 1) % size;

      rng ixr, oxr, yr(0, setup->ny - 1), zr(setup->nz - 1);

      ixr = rng( // input xr (halo)
        (i_min(setup->nx, rank, size) - halo + setup->nx) % setup->nx, 
        (i_min(setup->nx, rank, size) - 1    + setup->nx) % setup->nx  
      );
      oxr = rng( // output xr (edge)
        i_min(setup->nx, rank, size), 
        i_min(setup->nx, rank, size) + halo - 1
      );
      lhalo.reset(new slv_halo(this, peer_left, ixr, oxr, yr, zr));

      ixr = rng( // input xr (halo)
        (i_max(setup->nx, rank, size) + 1    + setup->nx) % setup->nx, 
        (i_max(setup->nx, rank, size) + halo + setup->nx) % setup->nx
      );
      oxr = rng( // output xr (edge)
        i_max(setup->nx, rank, size) - halo + 1,
        i_max(setup->nx, rank, size)  
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
