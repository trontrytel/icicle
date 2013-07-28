#include "kin_cloud_2d_common.hpp"

#include <libcloudph++/blk_2m/options.hpp>

template <
  typename real_t, 
  int n_iters, 
  typename ix,
  int n_eqs = 6
>
class kin_cloud_2d_blk_2m : public kin_cloud_2d_common<real_t, n_iters, ix, n_eqs>
{
  using parent_t = kin_cloud_2d_common<real_t, n_iters, ix, n_eqs>;

  //
  void update_forcings(libmpdataxx::arrvec_t<typename parent_t::arr_t> &rhs)
  {
    parent_t::update_forcings(rhs);
  }

  libcloudphxx::blk_2m::opts_t<real_t> opts;

  public:

  struct params_t : parent_t::params_t 
  { 
    libcloudphxx::blk_2m::opts_t<real_t> cloudph_opts;
  };

  // ctor
  kin_cloud_2d_blk_2m( 
    typename parent_t::ctor_args_t args, 
    const params_t &p
  ) : 
    parent_t(args, p),
    opts(p.cloudph_opts)
  { }  
};

