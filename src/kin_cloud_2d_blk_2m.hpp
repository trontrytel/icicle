#include "kin_cloud_2d_common.hpp"

#include <libcloudph++/blk_2m/options.hpp>
#include <libcloudph++/blk_2m/rhs_cellwise.hpp>
#include <libcloudph++/blk_2m/rhs_columnwise.hpp>

template <class ct_params_t>
class kin_cloud_2d_blk_2m : public kin_cloud_2d_common<ct_params_t>
{
  using parent_t = kin_cloud_2d_common<ct_params_t>;

  public:
  using ix = typename ct_params_t::ix;
  using real_t = typename ct_params_t::real_t;
  private:

  void zero_if_uninitialised(int e)  //TODO move to common
  {
    if (!finite(sum(this->state(e)(this->ijk)))) 
    this->state(e)(this->ijk) = 0;
  }

  void update_rhs(
    libmpdataxx::arrvec_t<typename parent_t::arr_t> &rhs,
    const typename parent_t::real_t &dt,
    const int &at 
  ) {
    parent_t::update_rhs(rhs, dt, at);

    this->mem->barrier(); // TODO: if neccesarry, then move to adv_rhs/....hpp

    // cell-wise
    {
      auto
	dot_rhod_th = rhs.at(ix::rhod_th)(this->i, this->j),
	dot_rhod_rv = rhs.at(ix::rhod_rv)(this->i, this->j),
	dot_rhod_rc = rhs.at(ix::rhod_rc)(this->i, this->j),
	dot_rhod_rr = rhs.at(ix::rhod_rr)(this->i, this->j),
	dot_rhod_nc = rhs.at(ix::rhod_nc)(this->i, this->j),
	dot_rhod_nr = rhs.at(ix::rhod_nr)(this->i, this->j);
      const auto
	rhod     = this->rhod(this->i, this->j),
	rhod_th  = this->state(ix::rhod_th)(this->i, this->j),
	rhod_rv  = this->state(ix::rhod_rv)(this->i, this->j),
	rhod_rc  = this->state(ix::rhod_rc)(this->i, this->j),
	rhod_rr  = this->state(ix::rhod_rr)(this->i, this->j),
	rhod_nc  = this->state(ix::rhod_nc)(this->i, this->j),
	rhod_nr  = this->state(ix::rhod_nr)(this->i, this->j);

      libcloudphxx::blk_2m::rhs_cellwise<real_t>(
        opts, dot_rhod_th, dot_rhod_rv, dot_rhod_rc, dot_rhod_nc, dot_rhod_rr, dot_rhod_nr,
	rhod,     rhod_th,     rhod_rv,     rhod_rc,     rhod_nc,     rhod_rr,     rhod_nr,
        this->dt
      );
    }

    // column-wise
    for (int i = this->i.first(); i <= this->i.last(); ++i)
    {
      auto
	dot_rhod_rr = rhs.at(ix::rhod_rr)(i, this->j),
	dot_rhod_nr = rhs.at(ix::rhod_nr)(i, this->j);
      const auto
        rhod     = this->rhod(i, this->j),
	rhod_rr  = this->state(ix::rhod_rr)(i, this->j),
	rhod_nr  = this->state(ix::rhod_nr)(i, this->j);

      libcloudphxx::blk_2m::rhs_columnwise<real_t>(
        opts, 
        dot_rhod_rr, dot_rhod_nr, 
        rhod,
            rhod_rr,     rhod_nr,  
	this->dt,
	this->dz
      );
    }

    this->mem->barrier(); // TODO: if needed, move to adv+rhs
  }

  libcloudphxx::blk_2m::opts_t<real_t> opts;

  protected:

  void hook_ante_loop(int nt)
  {
    // if uninitialised fill with zeros
    zero_if_uninitialised(ix::rhod_rc);
    zero_if_uninitialised(ix::rhod_nc);
    zero_if_uninitialised(ix::rhod_rr);
    zero_if_uninitialised(ix::rhod_nr);

    parent_t::hook_ante_loop(nt); // forcings after adjustments
  }

  // spinup stuff
  bool get_rain() { return opts.acnv; }
  void set_rain(bool val) 
  { 
    opts.acnv = val; 
    opts.RH_max = val ? 44 : 1.01; // 1% limit during spinup
  };

  public:

  struct rt_params_t : parent_t::rt_params_t 
  { 
    libcloudphxx::blk_2m::opts_t<real_t> cloudph_opts;
  };

  // ctor
  kin_cloud_2d_blk_2m( 
    typename parent_t::ctor_args_t args, 
    const rt_params_t &p
  ) : 
    parent_t(args, p),
    opts(p.cloudph_opts)
  { 
    assert(p.dt != 0);
  }  
};

