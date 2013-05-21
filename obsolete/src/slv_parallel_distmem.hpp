/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#pragma once
#include "slv_parallel.hpp" 

template <typename real_t, class shrdmem_class>
class slv_parallel_distmem : public shrdmem_class
{
  // <nested class>
  private: class slv_halo : public slv<real_t>
  {
    private: slv_parallel_distmem<real_t, shrdmem_class> *nghbr;
    private: int peer, cnt;
    private: unique_ptr<mtx::arr<real_t> > ibuf, obuf;
    private: mtx::rng ixr, oxr, yr, zr;

    public: slv_halo(slv_parallel_distmem<real_t, shrdmem_class> *nghbr, int peer,
      const mtx::rng &ixr, const mtx::rng &oxr, const mtx::rng &yr, const mtx::rng &zr
    )
      : nghbr(nghbr), peer(peer), ixr(ixr), oxr(oxr), yr(yr), zr(zr)
    { 
      ibuf.reset(new mtx::arr<real_t>(mtx::idx_ijk(ixr, yr, zr)));
      obuf.reset(new mtx::arr<real_t>(mtx::idx_ijk(oxr, yr, zr)));
      cnt = ibuf->cols() * ibuf->rows() * ibuf->depth();
    }

    public: typename mtx::arr<real_t>::type data(int e, int n, const mtx::idx &idx)
    { 
      return (*ibuf)(idx); 
    }
 
    public: void sync(int e, int n)
    {
      (*obuf)(ixr, yr, zr) = nghbr->data(e, n, mtx::idx_ijk(oxr, yr, zr));
      nghbr->sndrcv(peer, cnt, ibuf->data(), obuf->data());
    }

    public: void integ_loop() 
    { 
      assert(false); 
    }
  };
  // </nested class>

  private: int size, rank;
  private: unique_ptr<slv_halo> lhalo, rhalo;
  private: unique_ptr<mtx::arr<real_t> > libuf, ribuf, lobuf, robuf;

  private: int i_min(int nx, int rank, int size) { return (rank + 0) * nx / size; }
  private: int i_max(int nx, int rank, int size) { return i_min(nx, rank + 1, size) - 1; }

  // ctor
  public: slv_parallel_distmem(
    const stp<real_t> &setup, 
    out<real_t> &output,
    int size, 
    int rank
  )
    : shrdmem_class(setup, output,
      i_min(setup.grid.nx(), rank, size), i_max(setup.grid.nx(), rank, size), 
      0, setup.grid.ny() - 1,
      0, setup.grid.nz() - 1,
      1 // FIXME: MPI+OpenMP and MPI+threads 
    ), size(size), rank(rank)
  {
    // we distinguish left and right halos just by knowing where they come from so...
    if (size <= 2) error_macro("at least three subdomains are needed for distmem parallelisations")

    // subdomain length
    int nxs = setup.grid.nx() / size;
    if (nxs != ((1.*setup.grid.nx()) / (1.*size))) 
      error_macro("nx/nk must be an integer value (" << setup.grid.nx() << "/" << size << " given)")

    // halo containers for halo domain (TODO: use some clever Blitz constructor to save memory)
    {
      int halo = (setup.advsch.stencil_extent() - 1) / 2;
      int peer_left = (size + rank - 1) % size;
      int peer_rght = (size + rank + 1) % size;

      mtx::rng ixr, oxr, yr(0, setup.grid.ny() - 1), zr(setup.grid.nz() - 1);

      ixr = mtx::rng( // input xr (halo)
        (i_min(setup.grid.nx(), rank, size) - halo + setup.grid.nx()) % setup.grid.nx(), 
        (i_min(setup.grid.nx(), rank, size) - 1    + setup.grid.nx()) % setup.grid.nx()  
      );
      oxr = mtx::rng( // output xr (edge)
        i_min(setup.grid.nx(), rank, size), 
        i_min(setup.grid.nx(), rank, size) + halo - 1
      );
      lhalo.reset(new slv_halo(this, peer_left, ixr, oxr, yr, zr));

      ixr = mtx::rng( // input xr (halo)
        (i_max(setup.grid.nx(), rank, size) + 1    + setup.grid.nx()) % setup.grid.nx(), 
        (i_max(setup.grid.nx(), rank, size) + halo + setup.grid.nx()) % setup.grid.nx()
      );
      oxr = mtx::rng( // output xr (edge)
        i_max(setup.grid.nx(), rank, size) - halo + 1,
        i_max(setup.grid.nx(), rank, size)  
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
