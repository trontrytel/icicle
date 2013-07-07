#pragma once

#include "kin_cloud_2d_common.hpp"

#include <libcloudph++/lgrngn/options.hpp>


// TODO: what should the CPU do while waining for GPU working on particles stuff?

// @brief a minimalistic kinematic cloud model with lagrangian microphysics
//        built on top of the mpdata_2d solver (by extending it with
//        custom hook_ante_loop() and hook_post_step() methods)
template <
  typename real_t, 
  int n_iters, 
  solvers::inhomo_e inhomo,
  typename ix,
  int n_eqs = 4
>
class kin_cloud_2d_lgrngn : public kin_cloud_2d_common<real_t, n_iters, inhomo, ix, n_eqs>
{
  using parent_t = kin_cloud_2d_common<real_t, n_iters, inhomo, ix, n_eqs>; // TODO: lgrngn has no inhomo - just adjustments

  protected:

  // deals with initial supersaturation
  void hook_ante_loop(int nt)
  {
//    condevap(); // TODO: max supersaturation for spin-up?
    parent_t::hook_ante_loop(nt); // forcings after adjustments
  }

  // 
  void hook_post_step()
  {
// TODO: all the sdm logic
    parent_t::hook_post_step(); // includes the above forcings
  }

  real_t dx, dy, dz;
  libcloudphxx::lgrngn::opts<real_t> opts;

  public:

  struct params_t : parent_t::params_t 
  { 
    real_t dx = 0, dz = 0;
    libcloudphxx::lgrngn::opts<real_t> cloudph_opts;
  };

  // ctor
  kin_cloud_2d_lgrngn( 
    typename parent_t::ctor_args_t args, 
    const params_t &p
  ) : 
    parent_t(args, p),
    dx(p.dx),
    dy(1),
    dz(p.dz),
    opts(p.cloudph_opts)
  {
    assert(p.dx != 0);
    assert(p.dz != 0);
    assert(p.dt != 0); 
    opts.dt = p.dt;
  }  
};
