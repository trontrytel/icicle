#include "kin_cloud_2d_common.hpp"

#include <libcloudph++/blk_2m/options.hpp>
#include <libcloudph++/blk_2m/forcings_elementwise.hpp>

template <
  typename real_t, 
  int n_iters, 
  typename ix,
  int n_eqs = 6
>
class kin_cloud_2d_blk_2m : public kin_cloud_2d_common<real_t, n_iters, ix, n_eqs>
{
  using parent_t = kin_cloud_2d_common<real_t, n_iters, ix, n_eqs>;

  void zero_if_uninitialised(int e)  //TODO move to common
  {
    if (!finite(sum(this->state(e)(this->ijk)))) 
    this->state(e)(this->ijk) = 0;
  }

  void update_forcings(libmpdataxx::arrvec_t<typename parent_t::arr_t> &rhs)
  {
    parent_t::update_forcings(rhs);

    auto
      drhod_th = rhs.at(ix::rhod_th),       //TODO used just once -> probably an overkill
      drhod_rv = rhs.at(ix::rhod_rv),
      drhod_rc = rhs.at(ix::rhod_rc),
      drhod_rr = rhs.at(ix::rhod_rr),
      drhod_nc = rhs.at(ix::rhod_nc),
      drhod_nr = rhs.at(ix::rhod_nr);
    const auto
      rhod     = this->rhod,
      rhod_th  = this->state(ix::rhod_th),
      rhod_rv  = this->state(ix::rhod_rv),
      rhod_rc  = this->state(ix::rhod_rc),
      rhod_rr  = this->state(ix::rhod_rr),
      rhod_nc  = this->state(ix::rhod_nc),
      rhod_nr  = this->state(ix::rhod_nr);

    // element-wise
    {
      this->mem->barrier(); //TODO is it really needed?

      const rng_t &i = this->i, &j = this->j;
      libcloudphxx::blk_2m::forcings_elementwise<real_t>(opts,     drhod_th(i,j), drhod_rv(i,j), drhod_rc(i,j), drhod_nc(i,j), 
		                                         rhod(i,j), rhod_th(i,j),  rhod_rv(i,j),  rhod_rc(i,j),  rhod_nc(i,j));
      this->mem->barrier(); //TODO is it really needed?
    }

  }

  libcloudphxx::blk_2m::opts_t<real_t> opts;

  protected:

  void hook_ante_loop(int nt)
  {
    // if uninitialised fill with zeros
    zero_if_uninitialised(ix::rhod_rc);
    zero_if_uninitialised(ix::rhod_nc);
    zero_if_uninitialised(ix::rhod_rr);
    zero_if_uninitialised(ix::rhod_nc);

    parent_t::hook_ante_loop(nt); // forcings after adjustments
  }

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
  { 
    assert(p.dt != 0);
    opts.dt = p.dt;
  }  
};

