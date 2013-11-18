#include "kin_cloud_2d_common.hpp"

#include <libcloudph++/blk_1m/options.hpp>
#include <libcloudph++/blk_1m/adj_cellwise.hpp>
#include <libcloudph++/blk_1m/rhs_cellwise.hpp>
#include <libcloudph++/blk_1m/rhs_columnwise.hpp>

// @brief a minimalistic kinematic cloud model with bulk microphysics
//        built on top of the mpdata_2d solver (by extending it with
//        custom hook_ante_loop() and hook_post_step() methods)
template <
  typename real_t, 
  typename ix
>
class kin_cloud_2d_blk_1m : public kin_cloud_2d_common<real_t, ix>
{
  using parent_t = kin_cloud_2d_common<real_t, ix>;

  void condevap()
  {
    auto 
      rhod_th = this->state(ix::rhod_th)(this->ijk), // dry static energy density divided by c_pd (= dry air density times theta)
      rhod_rv = this->state(ix::rhod_rv)(this->ijk), // water vapour density
      rhod_rc = this->state(ix::rhod_rc)(this->ijk), // cloud water density
      rhod_rr = this->state(ix::rhod_rr)(this->ijk); // rain water density
    auto const
      rhod    = this->rhod(this->ijk);
      
    libcloudphxx::blk_1m::adj_cellwise<real_t>( 
      opts, rhod, rhod_th, rhod_rv, rhod_rc, rhod_rr, this->dt
    );
    this->mem->barrier(); 
  }

  void zero_if_uninitialised(int e)
  {
    if (!finite(sum(this->state(e)(this->ijk)))) 
      this->state(e)(this->ijk) = 0;
  }

  protected:

  bool get_rain() { return opts.conv; }
  void set_rain(bool val) { opts.conv = val; };

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

    // cell-wise
    {
      auto 
	dot_rhod_rc = rhs.at(ix::rhod_rc)(this->i, this->j),
	dot_rhod_rr = rhs.at(ix::rhod_rr)(this->i, this->j);
      const auto 
	rhod_rc  = this->state(ix::rhod_rc)(this->i, this->j),
	rhod_rr  = this->state(ix::rhod_rr)(this->i, this->j),
	rhod     = this->rhod(this->i, this->j);
      libcloudphxx::blk_1m::rhs_cellwise<real_t>(opts, dot_rhod_rc, dot_rhod_rr, rhod, rhod_rc, rhod_rr);
    }

    // column-wise
    for (int i = this->i.first(); i <= this->i.last(); ++i)
    {
      auto 
	dot_rhod_rr = rhs.at(ix::rhod_rr)(i, this->j);
      const auto 
	rhod_rr  = this->state(ix::rhod_rr)(i, this->j),
	rhod     = this->rhod(i, this->j);
      libcloudphxx::blk_1m::rhs_columnwise<real_t>(opts, dot_rhod_rr, rhod, rhod_rr, this->dz);
    }
  }

  // 
  void hook_post_step()
  {
    condevap(); // treat saturation adjustment as post-advection, pre-rhs adjustment
    parent_t::hook_post_step(); // includes the above forcings
  }

  libcloudphxx::blk_1m::opts_t<real_t> opts;

  public:

  struct params_t : parent_t::params_t 
  { 
    libcloudphxx::blk_1m::opts_t<real_t> cloudph_opts;
  
    // ctor
    params_t() { this->n_eqs = 4; }
  };

  // ctor
  kin_cloud_2d_blk_1m( 
    typename parent_t::ctor_args_t args, 
    const params_t &p
  ) : 
    parent_t(args, p),
    opts(p.cloudph_opts)
  {}  
};
