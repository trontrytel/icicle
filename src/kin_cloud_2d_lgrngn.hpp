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
  libcloudphxx::lgrngn::opts_t<real_t> opts;
  std::unique_ptr<libcloudphxx::lgrngn::particles_proto<real_t>> prtcls;

  void diag()
  {
    assert(this->mem->rank() == 0);
    prtcls->diag();
    // TODO: specify somehow what to record... and pass it to record_aux()
    this->record_aux("sd_conc", prtcls->outbuf());
    //                          ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ TODO!!!
  } 

  protected:

  // deals with initial supersaturation
  void hook_ante_loop(int nt)
  {
    // TODO: max supersaturation for spin-up?
    parent_t::hook_ante_loop(nt); 

    // TODO: barrier?
    if (this->mem->rank() == 0) 
    {
      this->setup_aux("sd_conc"); // TODO: setup other fields, units?

      prtcls->init(
        this->rhod.stride().data(),
	this->mem->state(ix::rhod_th).dataZero(),
	this->mem->state(ix::rhod_rv).dataZero(),
	this->rhod.dataZero()
      ); 
      diag();
    }
    // TODO: barrier?
  }

  // 
  void hook_post_step()
  {
    parent_t::hook_post_step(); // includes output
    // TODO: barrier?
    if (this->mem->rank() == 0) 
    {
      prtcls->step(
	this->state(ix::rhod_th).dataZero(),
	this->state(ix::rhod_rv).dataZero()
      ); 
      if (this->timestep % this->outfreq == 0) diag();
    }
    // TODO: barrier?
  }

  public:

  struct params_t : parent_t::params_t 
  { 
    int backend = -1;
    libcloudphxx::lgrngn::opts_t<real_t> cloudph_opts;
  };

  // ctor
  kin_cloud_2d_lgrngn( 
    typename parent_t::ctor_args_t args, 
    const params_t &p
  ) : 
    parent_t(args, p),
    opts(p.cloudph_opts)
  {
    assert(p.backend != -1);

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
