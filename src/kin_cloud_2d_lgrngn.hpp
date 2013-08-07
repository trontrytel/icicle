#pragma once

#include "kin_cloud_2d_common.hpp"
#include "outmom.hpp"

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
  using parent_t = kin_cloud_2d_common<real_t, n_iters, ix, n_eqs>; // TODO: lgrngn has no inhomo - just adjustments

  // member fields
  std::unique_ptr<libcloudphxx::lgrngn::particles_proto<real_t>> prtcls;
  outmom_t<real_t> out_dry, out_wet;

  // helper methods
  void diag()
  {
    assert(this->mem->rank() == 0);

    // recording super-droplet concentration per grid cell 
    prtcls->diag_sd_conc();
    this->record_aux("sd_conc", prtcls->outbuf());
   
    // recording requested statistical moments
    {
      // dry
      int rng_num = 0;
      for (auto &rng_moms : out_dry)
      {
        auto &rng(rng_moms.first);
        prtcls->diag_dry_rng(rng.first / si::metres, rng.second / si::metres);
        for (auto &mom : rng_moms.second)
        {
          prtcls->diag_dry_mom(mom);
          this->record_aux(aux_name("rd", rng_num, mom), prtcls->outbuf());
        }
        rng_num++;
      }
    }
    {
      // wet
      int rng_num = 0;
      for (auto &rng_moms : out_wet)
      {
        auto &rng(rng_moms.first);
        prtcls->diag_wet_rng(rng.first / si::metres, rng.second / si::metres);
        for (auto &mom : rng_moms.second)
        {
          prtcls->diag_wet_mom(mom);
          this->record_aux(aux_name("rw", rng_num, mom), prtcls->outbuf());
        }
        rng_num++;
      }
    }
  } 

  libcloudphxx::lgrngn::arrinfo_t<real_t> make_arrinfo(
    typename parent_t::arr_t arr
  ) {
    return libcloudphxx::lgrngn::arrinfo_t<real_t>(
      arr.dataZero(), 
      arr.stride().data()
    );
  }

  std::string aux_name(
    const std::string pfx, 
    const int rng,
    const int mom
  )
  { 
    std::ostringstream tmp;
    tmp << pfx << "_rng" << rng << "_mom" << mom;
    return tmp.str();
  }

  void setup_aux_helper(
    const std::string pfx, 
    const outmom_t<real_t> &outmoms
  )
  {
    // TODO: attributes incl. units?
    int rng = 0; // TODO: perhaps separate numbering for rd and rw - not to mislead?
    for (auto &rng_moms : outmoms)
    {
      for (auto &mom : rng_moms.second)
      {
	this->setup_aux(aux_name(pfx, rng, mom)); 
      }
      ++rng;
    }
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
      setup_aux_helper("rd", out_dry); 
      setup_aux_helper("rw", out_wet); 
      this->setup_aux("sd_conc"); 

      prtcls->init(
        make_arrinfo(this->mem->state(ix::rhod_th)),
        make_arrinfo(this->mem->state(ix::rhod_rv)),
	make_arrinfo(this->rhod),
        make_arrinfo(this->mem->courant(0)),
        make_arrinfo(this->mem->courant(1))
      ); 

      // writing diagnostic data for the initial condition
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
    outmom_t<real_t> out_dry, out_wet;
  };

  // ctor
  kin_cloud_2d_lgrngn( 
    typename parent_t::ctor_args_t args, 
    const params_t &p
  ) : 
    parent_t(args, p),
    out_dry(p.out_dry), // TODO: unnecesarily per each thread!
    out_wet(p.out_wet)  // ditto
  {
    if (this->mem->rank() == 0)
    {
      assert(p.backend != -1);
      assert(p.dt != 0); 

      auto opts(p.cloudph_opts); // making a copy to be able to modify it
      opts.dt = p.dt; // advection timestep = microphysics timestep
      opts.dx = p.dx;
      opts.dz = p.dz;

      prtcls.reset(libcloudphxx::lgrngn::factory<real_t>::make(p.backend, opts));
    }
  }  
};
