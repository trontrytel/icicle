#pragma once

#include "kin_cloud_2d_common.hpp"
#include "outmom.hpp"

#include <libcloudph++/lgrngn/options.hpp>
#include <libcloudph++/lgrngn/particles.hpp>


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
      for (auto &rng_moms : params.out_dry)
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
      for (auto &rng_moms : params.out_wet)
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
    parent_t::hook_ante_loop(nt); 

    // TODO: barrier?
    if (this->mem->rank() == 0) 
    {
      assert(params.backend != -1);
      assert(params.dt != 0); 

      params.cloudph_opts.dt = params.dt; // advection timestep = microphysics timestep
      params.cloudph_opts.dx = params.dx;
      params.cloudph_opts.dz = params.dz;

      prtcls.reset(libcloudphxx::lgrngn::factory<real_t>::make(params.backend, params.cloudph_opts));

      setup_aux_helper("rd", params.out_dry); 
      setup_aux_helper("rw", params.out_wet); 
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
      // assuring previous async step finished (but not in first timestep)
      if (params.async && ftr.valid()) ftr.get();

      // running synchronous stuff
      prtcls->step_sync(
        make_arrinfo(this->mem->state(ix::rhod_th)),
        make_arrinfo(this->mem->state(ix::rhod_rv))
      ); 

      // running asynchronous stuff
      {
        using libcloudphxx::lgrngn::particles;
        using libcloudphxx::lgrngn::cuda;

        if (params.async && params.backend == cuda)
        {
          ftr = std::async(
            std::launch::async, 
            &particles<real_t, cuda>::step_async, 
            dynamic_cast<particles<real_t, cuda>*>(prtcls.get())
          );
        } else prtcls->step_async();
      }

      // performing diagnostics
      if (this->timestep % this->outfreq == 0) 
      { 
        if (params.async)
        {
          assert(ftr.valid());
          ftr.get();
        }
        diag();
      }
    }
    // TODO: barrier?
  }

  public:

  struct params_t : parent_t::params_t 
  { 
    int backend = -1;
    bool async = true;
    libcloudphxx::lgrngn::opts_t<real_t> cloudph_opts;
    outmom_t<real_t> out_dry, out_wet;
  };

  private:

  // per-thread copy of params
  params_t params;

  public:

  // ctor
  kin_cloud_2d_lgrngn( 
    typename parent_t::ctor_args_t args, 
    const params_t &p
  ) : 
    parent_t(args, p),
    params(p)
  {
    // delaying any initialisation to ante_loop as rank() does not function within ctor!
    // TODO: equip rank() in libmpdata with an assert() checking if not in serial block
  }  
};
