#pragma once

#include "kin_cloud_2d_common.hpp"

#include <libcloudph++/lgrngn/options.hpp>
#include <libcloudph++/lgrngn/particles.hpp>


// TODO: what should the CPU do while waining for GPU working on particles stuff?

// @brief a minimalistic kinematic cloud model with lagrangian microphysics
//        built on top of the mpdata_2d solver (by extending it with
//        custom hook_ante_loop() and hook_post_step() methods)
template <
  typename real_t, 
  int n_iters, 
  typename ix,
  int n_eqs = 4
>
class kin_cloud_2d_lgrngn : public kin_cloud_2d_common<real_t, n_iters, ix, n_eqs>
{
  using parent_t = kin_cloud_2d_common<real_t, n_iters, ix, n_eqs>; 
  // TODO: lgrngn has no inhomo - just adjustments

  // member fields
  real_t dx, dz;
  libcloudphxx::lgrngn::opts_t<real_t> opts;
  std::unique_ptr<libcloudphxx::lgrngn::particles_proto<real_t>> prtcls;

  protected:

  // deals with initial supersaturation
  void hook_ante_loop(int nt)
  {
    // TODO: max supersaturation for spin-up?
    parent_t::hook_ante_loop(nt); // forcings after adjustments
    // TODO: barrier?
    if (this->mem->rank() == 0) 
    {
      prtcls->init(
        this->rhod.stride().data(),
	this->mem->state(ix::rhod_th).dataZero(),
	this->mem->state(ix::rhod_rv).dataZero(),
	this->rhod.dataZero()
      ); 
    }
    // TODO: barrier?
  }

  // 
  void hook_post_step()
  {
    parent_t::hook_post_step();
    // TODO: barrier?
    if (this->mem->rank() == 0) 
    {
      prtcls->step(
	this->state(ix::rhod_th).dataZero(),
	this->state(ix::rhod_rv).dataZero()
      ); 
    }
    // TODO: barrier?
  }

  public:

  struct params_t : parent_t::params_t 
  { 
    real_t dx = 0, dz = 0;
    int backend = -1;
    libcloudphxx::lgrngn::opts_t<real_t> cloudph_opts;
  };

  // ctor
  kin_cloud_2d_lgrngn( 
    typename parent_t::ctor_args_t args, 
    const params_t &p
  ) : 
    parent_t(args, p),
    dx(p.dx),
    dz(p.dz),
    opts(p.cloudph_opts)
  {
    assert(p.dx != 0);
    assert(p.dz != 0);
    assert(p.dt != 0); 
    assert(p.backend != -1);
    //opts.dt = p.dt;

    prtcls.reset(libcloudphxx::lgrngn::factory<real_t>::make(
      p.backend, 
      p.cloudph_opts.sd_conc_mean,
      p.cloudph_opts.dry_distros,
      p.cloudph_opts.nx,
      p.dx,
      p.cloudph_opts.nz,
      p.dz      
    ));
  }  
};
