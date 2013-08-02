#pragma once

// TODO: other includes?
#include <libcloudph++/common/hydrostatic.hpp>
#include <libcloudph++/common/theta_std.hpp>
#include <libcloudph++/common/lognormal.hpp>
#include <libcloudph++/common/unary_function.hpp>

// TODO: relaxation terms still missing

// 8th ICMW case 1 by Wojciech Grabowski)
namespace icmw8_case1
{
  using real_t = double;

  namespace hydrostatic = libcloudphxx::common::hydrostatic;
  namespace theta = libcloudphxx::common::theta_std;
  namespace lognormal = libcloudphxx::common::lognormal;

  enum {x, z}; // dimensions

  const quantity<si::temperature, real_t> 
    th_0 = 289 * si::kelvins;
  const quantity<si::dimensionless, real_t> 
    rv_0 = 7.5e-3;
  const quantity<si::pressure, real_t> 
    p_0 = 101500 * si::pascals;
  const quantity<si::velocity, real_t> 
    w_max = real_t(.6) * si::metres_per_second;
  const quantity<si::length, real_t> 
    z_0  = 0    * si::metres,
    nzdz = 1500 * si::metres, 
    nxdx = 1500 * si::metres;
  const quantity<si::time, real_t>
    dt = .1 * si::seconds;//4 * si::seconds;

  //aerosol bimodal lognormal dist. for super droplets
  const quantity<si::length, real_t>
    mean_rd1 = .04e-6 / 2 * si::metres,
    mean_rd2 = .15e-6 / 2 * si::metres;
  const quantity<si::dimensionless, real_t>
    sdev_rd1 = 1.4,
    sdev_rd2 = 1.6;
  const quantity<power_typeof_helper<si::length, static_rational<-3>>::type, real_t>
    n1_tot = 60e6 / si::cubic_metres,
    n2_tot = 40e6 / si::cubic_metres;

  //aerosol lognormal dist. for 2-moment bulk scheme (needed for activation)
  auto mean_rd = mean_rd2;
  auto sdev_rd = sdev_rd2;
  auto N_tot   = n2_tot;
  //aerosol chemical composition parameter (needed for activation)
  const quantity<si::dimensionless, real_t> chem_b = .55; //ammonium sulphate
                                          //chem_b = 1.33; // sodium chloride

  // density profile as a function of altitude
  struct rhod
  {
    real_t operator()(real_t z) const
    {
      quantity<si::pressure, real_t> p = hydrostatic::p(
	z * si::metres, th_0, rv_0, z_0, p_0
      );
      
      quantity<si::mass_density, real_t> rhod = theta::rhod(
	p, th_0, rv_0
      );

      return rhod / si::kilograms * si::cubic_metres;
    }

    // to make the rhod() functor accept Blitz arrays as arguments
    BZ_DECLARE_FUNCTOR(rhod);
  };



  /// (similar to eq. 2 in @copydetails Rasinski_et_al_2011, Atmos. Res. 102)
  /// @arg xX = x / (nx*dx)
  /// @arg zZ = z / (nz*dz)
  real_t psi(real_t xX, real_t zZ)
  {
    return - sin(pi<real_t>() * zZ) * cos(2 * pi<real_t>() * xX);
  }
  BZ_DECLARE_FUNCTION2_RET(psi, real_t)

  // function expecting a libmpdata solver parameters struct as argument
  template <class T>
  void setopts(T &params, int nx, int nz)
  {
    params.dt = dt / si::seconds;
    params.dx = nxdx / si::metres / nx;
    params.dz = nzdz / si::metres / nz;

    params.rhod.resize(nz);
    for (int j = 0; j < nz; ++j)
      params.rhod[j] = rhod()((j+.5) * params.dz);
  }

  // function expecting a libmpdata++ solver as argument
  template <class concurr_t>
  void intcond(concurr_t &solver)
  {
    using ix = typename concurr_t::solver_t::ix;

    // helper ondex placeholders
    blitz::firstIndex i;
    blitz::secondIndex j;

    // dx, dy ensuring 1500x1500 domain
    int 
      nx = solver.state().extent(x),
      nz = solver.state().extent(z);
    real_t 
      dx = nxdx / si::metres / nx, 
      dz = nzdz / si::metres / nz; 

    // constant potential temperature & water vapour mixing ratio profiles
    solver.state(ix::rhod_th) = rhod()((j+.5)*dz) * (th_0 / si::kelvins); // TODO: should be theta_dry and is theta
    solver.state(ix::rhod_rv) = rhod()((j+.5)*dz) * real_t(rv_0);

    // velocity field obtained by numerically differentiating a stream function
    {
      real_t A = (w_max / si::metres_per_second) * solver.state().extent(x) * dx / pi<real_t>();

      solver.courant(x) = - A * (
	psi(i/real_t(nx), (j+.5+.5)/nz)-
	psi(i/real_t(nx), (j+.5-.5)/nz)
      ) / dz             // numerical derivative
      / rhod()((j+.5)* dz) // psi defines rho_d times velocity
      * (dt / si::seconds) / dx;         // converting to Courant number

      solver.courant(z) = A * (
	psi((i+.5+.5)/nx, j/real_t(nz)) -
	psi((i+.5-.5)/nx, j/real_t(nz))
      ) / dx 
      / rhod()(j * dz)
      * (dt / si::seconds) / dz; 
    }
  }

  // lognormal aerosol distribution
  template <typename T>
  struct log_dry_radii : public libcloudphxx::common::unary_function<T>
  {
    T funval(const T lnrd) const
    {
      return T((
          lognormal::n_e(mean_rd1, sdev_rd1, n1_tot, quantity<si::dimensionless, real_t>(lnrd)) +
          lognormal::n_e(mean_rd2, sdev_rd2, n2_tot, quantity<si::dimensionless, real_t>(lnrd)) 
        ) * si::cubic_metres
      );
    }

    log_dry_radii *do_clone() const 
    { return new log_dry_radii( *this ); }
  };
};
