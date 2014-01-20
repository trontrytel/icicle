#pragma once

#include <libmpdata++/solvers/mpdata_rhs.hpp>
#include <libmpdata++/output/hdf5.hpp>

using namespace libmpdataxx; // TODO: get rid of it?

template <class ct_params_t>
class kin_cloud_2d_common : public 
  output::hdf5<
    solvers::mpdata_rhs<ct_params_t>
  >
{
  using parent_t = output::hdf5<solvers::mpdata_rhs<ct_params_t>>;

  protected:

  typename parent_t::arr_t rhod;
  typename ct_params_t::real_t dx, dz; // 0->dx, 1->dy ! TODO
  int spinup; // number of timesteps

  // spinup stuff
  virtual bool get_rain() = 0;
  virtual void set_rain(bool) = 0;

  void hook_ante_loop(int nt) 
  {
    if (get_rain() == false) spinup = 0; // spinup does not make sense without autoconversion  (TODO: issue a warning?)
    if (spinup > 0) set_rain(false);

    parent_t::hook_ante_loop(nt); 
  }

  void hook_ante_step()
  {
    // turn autoconversion on only after spinup (if spinup was specified)
    if (spinup != 0 && spinup == this->timestep) set_rain(true);

    parent_t::hook_ante_step(); 
  }


  void update_rhs(
    arrvec_t<typename parent_t::arr_t> &rhs,
    const typename parent_t::real_t &dt,
    const int &at 
  )   
  {   
    parent_t::update_rhs(rhs, dt, at);

    // relaxation terms
    {
      // computed level-wise
      //for (int j = 0; j < this->span[0]; ++j)
      {  
        // rhod_th
        {
          // TODO ...
          //const auto &psi = this->state(ix::rhod_th);
          //const auto psi_mean = this->mem->sum(psi, this->i, j) / this->span[0];
          //rhs.at(ix::rhod_th)(this->i, j) = - (psi(this->i, j) - psi_mean) / tau(j...);
        }
        // rhod_rv
        {
          // TODO: same as above...
        }
      }
    }

  }

  public:

  struct rt_params_t : parent_t::rt_params_t 
  { 
    std::vector<typename ct_params_t::real_t> rhod; // profile
    typename ct_params_t::real_t dx = 0, dz = 0;
    int spinup = 0; // number of timesteps during which autoconversion is to be turned off
  };

  // ctor
  kin_cloud_2d_common( 
    typename parent_t::ctor_args_t args, 
    const rt_params_t &p
  ) : 
    parent_t(args, p),
    rhod(args.mem->tmp[__FILE__][0][0]),
    dx(p.dx),
    dz(p.dz),
    spinup(p.spinup)
  {
    assert(p.rhod.size() == this->j.last()+1);
    assert(dx != 0);
    assert(dz != 0);

    // initialising rhod array columnwise with data from the p.rhod profile
    for (int i = this->i.first(); i <= this->i.last(); ++i)
      for (int j = this->j.first(); j <= this->j.last(); ++j)
	rhod(i, j) = p.rhod[j];
  }  

  static void alloc(typename parent_t::mem_t *mem, const rt_params_t &p)
  {
    using namespace libmpdataxx::arakawa_c;
    parent_t::alloc(mem, p);
    const rng_t i(0, p.span[0]-1), j(0, p.span[1]-1);
    mem->tmp[__FILE__].push_back(new arrvec_t<typename parent_t::arr_t>());
    mem->tmp[__FILE__].back().push_back(new typename parent_t::arr_t(i^parent_t::halo, j^parent_t::halo)); // rhod
  }
};
