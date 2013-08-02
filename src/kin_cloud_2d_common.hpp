#pragma once

#include <libmpdata++/solvers/adv/mpdata_2d.hpp>
// note: FCT cannot be used as of now as the density is not constant in space here!
#include <libmpdata++/solvers/adv/mpdata_fct_2d.hpp>
#include <libmpdata++/solvers/adv+rhs/solver_inhomo.hpp>
#include <libmpdata++/output/hdf5.hpp>

using namespace libmpdataxx; // TODO: not here?

template <
  typename real_t,
  int n_iters,
  typename ix_t,
  int n_eqs
>
class kin_cloud_2d_common : public 
    output::hdf5<
      solvers::inhomo_solver<
	solvers::mpdata_2d<real_t, n_iters, n_eqs>, 
	solvers::strang
    >
  >
{
  using parent_t = output::hdf5<solvers::inhomo_solver<solvers::mpdata_2d<real_t, n_iters, n_eqs>, solvers::strang>>;

  protected:

  typename parent_t::arr_t rhod;
  real_t dx, dz; // 0->dx, 1->dy ! TODO

  public:

  typedef ix_t ix;

  struct params_t : parent_t::params_t 
  { 
    std::vector<real_t> rhod; // profile
    real_t dx = 0, dz = 0;
  };

  // ctor
  kin_cloud_2d_common( 
    typename parent_t::ctor_args_t args, 
    const params_t &p
  ) : 
    parent_t(args, p),
    rhod(args.mem->tmp[__FILE__][0][0]),
    dx(p.dx),
    dz(p.dz)
  {
    assert(p.rhod.size() == this->j.last()+1);
    assert(dx != 0);
    assert(dz != 0);

    // initialising rhod array columnwise with data from the p.rhod profile
    for (int i = this->i.first(); i <= this->i.last(); ++i)
      for (int j = this->j.first(); j <= this->j.last(); ++j)
	rhod(i, j) = p.rhod[j];
  }  

  static void alloc(typename parent_t::mem_t *mem, const int nx, const int ny)
  {
    using namespace libmpdataxx::arakawa_c;
    parent_t::alloc(mem, nx, ny);
    const rng_t i(0, nx-1), j(0, ny-1);
    mem->tmp[__FILE__].push_back(new arrvec_t<typename parent_t::arr_t>());
    mem->tmp[__FILE__].back().push_back(new typename parent_t::arr_t(i^parent_t::halo, j^parent_t::halo)); // rhod
  }
};
