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
  std::unique_ptr<libcloudphxx::lgrngn::particles_proto<real_t>> prtcls;

  void diag()
  {
    assert(this->mem->rank() == 0);
    prtcls->diag();
    // TODO: specify somehow what to record... and pass it to record_aux()
    this->record_aux("sd_conc", prtcls->outbuf());
    //                          ^^^^^^^^^^^^^^^^ TODO!!!
  } 

  libcloudphxx::lgrngn::arrinfo_t<real_t> make_arrinfo(
    typename parent_t::arr_t arr
  ) {
    return libcloudphxx::lgrngn::arrinfo_t<real_t>(
      arr.dataZero(), 
      arr.stride().data()
    );
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
        make_arrinfo(this->mem->state(ix::rhod_th)),
        make_arrinfo(this->mem->state(ix::rhod_rv)),
	make_arrinfo(this->rhod),
        make_arrinfo(this->mem->courant(0)),
        make_arrinfo(this->mem->courant(1))
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
        make_arrinfo(this->mem->state(ix::rhod_th)),
        make_arrinfo(this->mem->state(ix::rhod_rv))
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
    parent_t(args, p)
  {
    assert(p.backend != -1);

// TODO: repeated elswhere, and works on a copy... 
    assert(p.dt != 0); 
    auto opts = p.cloudph_opts;
    opts.dt = p.dt; // advection timestep = microphysics timestep

    prtcls.reset(libcloudphxx::lgrngn::factory<real_t>::make(p.backend, opts));
  }  
};
