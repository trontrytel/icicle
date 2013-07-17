#include "kin_cloud_2d_common.hpp"

#include <libcloudph++/blk_2m/options.hpp>

template <
  typename real_t, 
  int n_iters, 
  solvers::inhomo_e inhomo,
  typename ix,
  int n_eqs = 6
>
class kin_cloud_2d_blk_2m : public kin_cloud_2d_common<real_t, n_iters, inhomo, ix, n_eqs>
{
  using parent_t = kin_cloud_2d_common<real_t, n_iters, inhomo, ix, n_eqs>;

  //
  void update_forcings(libmpdataxx::arrvec_t<typename parent_t::arr_t> &rhs)
  {
    parent_t::update_forcings(rhs);
  }

  real_t dz; // TODO: all the dz stuff could perhaps go to common?
  libcloudphxx::blk_2m::opts<real_t> opts;

  public:

  struct params_t : parent_t::params_t 
  { 
    real_t dx = 0, dz = 0; // TODO: move to common?
    libcloudphxx::blk_2m::opts<real_t> cloudph_opts;
  };

  // ctor
  kin_cloud_2d_blk_2m( 
    typename parent_t::ctor_args_t args, 
    const params_t &p
  ) : 
    parent_t(args, p),
    dz(p.dz),
    opts(p.cloudph_opts)
  { 
    assert(p.dz != 0);
  }  
};

