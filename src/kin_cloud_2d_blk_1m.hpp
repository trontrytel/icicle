#include "kin_cloud_2d_common.hpp"

#include <libcloudph++/blk_1m/options.hpp>
#include <libcloudph++/blk_1m/adjustments.hpp>
#include <libcloudph++/blk_1m/forcings_elementwise.hpp>
#include <libcloudph++/blk_1m/forcings_columnwise.hpp>

// @brief a minimalistic kinematic cloud model with bulk microphysics
//        built on top of the mpdata_2d solver (by extending it with
//        custom hook_ante_loop() and hook_post_step() methods)
template <
  typename real_t, 
  int n_iters, 
  solvers::inhomo_e inhomo,
  typename ix,
  int n_eqs = 4
>
class kin_cloud_2d_blk_1m : public kin_cloud_2d_common<real_t, n_iters, inhomo, ix, n_eqs>
{
  using parent_t = kin_cloud_2d_common<real_t, n_iters, inhomo, ix, n_eqs>;

  void condevap()
  {
    auto 
      rhod_th = this->state(ix::rhod_th)(this->ijk), // dry static energy density divided by c_pd (= dry air density times theta)
      rhod_rv = this->state(ix::rhod_rv)(this->ijk), // water vapour density
      rhod_rc = this->state(ix::rhod_rc)(this->ijk), // cloud water density
      rhod_rr = this->state(ix::rhod_rr)(this->ijk); // rain water density
    auto const
      rhod    = this->rhod(this->ijk);
      
    libcloudphxx::blk_1m::adjustments<real_t>( 
      opts, rhod, rhod_th, rhod_rv, rhod_rc, rhod_rr
    );
  }

  void zero_if_uninitialised(int e)
  {
    if (!finite(sum(this->state(e)(this->ijk)))) 
      this->state(e)(this->ijk) = 0;
  }

  protected:

  // deals with initial supersaturation
  void hook_ante_loop(int nt)
  {
    // if uninitialised fill with zeros
    zero_if_uninitialised(ix::rhod_rc);
    zero_if_uninitialised(ix::rhod_rr);

    // deal with initial supersaturation
    condevap();

    parent_t::hook_ante_loop(nt); // forcings after adjustments
  }

  //
  void update_forcings(libmpdataxx::arrvec_t<typename parent_t::arr_t> &rhs)
  {
    parent_t::update_forcings(rhs);
 
    auto 
      drhod_rc = rhs.at(ix::rhod_rc),
      drhod_rr = rhs.at(ix::rhod_rr);
    const auto 
      rhod_rc  = this->state(ix::rhod_rc),
      rhod_rr  = this->state(ix::rhod_rr),
      rhod     = this->rhod;

    // element-wise
    {
      const rng_t &i = this->i, &j = this->j;
      libcloudphxx::blk_1m::forcings_elementwise<real_t>(opts, drhod_rc(i,j), drhod_rr(i,j), rhod(i,j), rhod_rc(i,j), rhod_rr(i,j));
    }

    // column-wise
    {
      const rng_t j = this->j;
      for (int i = this->i.first(); i <= this->i.last(); ++i)
	libcloudphxx::blk_1m::forcings_columnwise<real_t>(opts, drhod_rr(i,j), rhod(i,j), rhod_rr(i,j), dz);
    }
  }

  // 
  void hook_post_step()
  {
    condevap(); // treat saturation adjustment as post-advection, pre-rhs adjustment
    parent_t::hook_post_step(); // includes the above forcings
    // TODO: shouldn't condevap() be called again here to ensure adjusted field is output?
  }

  real_t dz;
  libcloudphxx::blk_1m::opts<real_t> opts;

  public:

  struct params_t : parent_t::params_t 
  { 
    real_t dz = 0;
    libcloudphxx::blk_1m::opts<real_t> cloudph_opts;
  };

  // ctor
  kin_cloud_2d_blk_1m( 
    typename parent_t::ctor_args_t args, 
    const params_t &p
  ) : 
    parent_t(args, p),
    dz(p.dz),
    opts(p.cloudph_opts)
  {
    assert(p.dz != 0);
    assert(p.dt != 0);
    opts.dt = p.dt;
  }  
};
