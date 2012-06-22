/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date March-June 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#pragma once

#include "sdm_ode.hpp"
#include "phc_theta.hpp"
#include "phc_const_cp.hpp"
#include "phc_kappa_koehler.hpp"

namespace sdm
{
  // identity: xi = rw
  template <typename T> struct xi_id 
  { 
    static T xi_of_rw(const T &rw) { return rw; }
    static T xi_of_rw3(const T &rw3) { return pow(rw3, T(1./3)); }
    static T rw_of_xi(const T &xi) { return xi; } 
    static T rw2_of_xi(const T &xi) { return xi * xi; } 
    static T rw3_of_xi(const T &xi) { return xi * xi * xi; } 
    static T dxidrw(const T &) { return 1; }
  };

  // xi = ln(rw / nm)
  template <typename T> struct xi_ln
  { 
    static T xi_of_rw(const T &rw) { return log(rw / T(1e-9)); }
    static T xi_of_rw3(const T &rw3) { return log(pow(rw3, T(1./3)) / T(1e-9)); }
    static T rw_of_xi(const T &xi) { return T(1e-9) * exp(xi); } 
    static T rw2_of_xi(const T &xi) { return pow(T(1e-9) * exp(xi), 2); } 
    static T rw3_of_xi(const T &xi) { return pow(T(1e-9) * exp(xi), 3); } 
    static T dxidrw(const T &rw) { return T(1) / rw; }
  };

  // xi = pow(rw, 2) 
  template <typename T> struct xi_p2
  { 
    static T xi_of_rw(const T &rw) { return rw * rw; }
    static T xi_of_rw3(const T &rw3) { return pow(rw3, T(2./3)); }
    static T rw_of_xi(const T &xi) { return sqrt(xi); } 
    static T rw2_of_xi(const T &xi) { return xi; } 
    static T rw3_of_xi(const T &xi) { return pow(xi, T(3./2)); } 
    static T dxidrw(const T &rw) { return T(2) * rw; }
  };

  // xi = pow(rw, 3) 
  template <typename T> struct xi_p3
  { 
    static T xi_of_rw(const T &rw) { return rw * rw * rw; }
    static T xi_of_rw3(const T &rw3) { return rw3; }
    static T rw_of_xi(const T &xi) { return pow(xi, T(1)/T(3)); } 
    static T rw2_of_xi(const T &xi) { return pow(xi, 2./3); } 
    static T rw3_of_xi(const T &xi) { return xi; } 
    static T dxidrw(const T &rw) { return T(3) * rw * rw; }
  };

  // nested the ODE system to be solved to update super-droplet sizes
  template <typename real_t, class algo, class xi>
  class ode_xi : public ode_algo<real_t, algo>
  { 
/*
    public: real_t transform(const real_t &x) const
    {   
      return xi::xi_of_rw(x);
    }  
*/

    // nested functor: setting wet radii to equilibrium values
    struct equil : xi
    { 
      stat_t<real_t> &stat;
      const envi_t<real_t> &envi;

      // ctor
      equil(
        stat_t<real_t> &stat, 
        const envi_t<real_t> &envi
      ) : stat(stat), envi(envi)
      { }

      void operator()(thrust_size_t id) 
      { 
        thrust_size_t ij = stat.ij[id]; 
        stat.xi[id] = this->xi_of_rw3(phc::rw3_eq<real_t>( // TODO: allow choice among different Koehler curves
          stat.rd3[id] * si::cubic_metres, 
          0 + stat.kpa[id], // it fails to compile without the zero!
          thrust::min(
            real_t(.99),
            real_t( // TODO: interpolation to drop positions?
              phc::R_v<real_t>() * (envi.T[ij] * si::kelvins) * (envi.rhod_rv[ij] * si::kilograms / si::cubic_metres) 
              / (phc::p_vs<real_t>(envi.T[ij] * si::kelvins))
            )
          )
        ) / si::cubic_metres); 
      }
    };

    // nested functor: the drop growth law
    private: 
    class drop_growth_equation : public xi
    {
      // private fields
      private: const stat_t<real_t> &stat;
      private: const envi_t<real_t> &envi;

      // ctor
      public: drop_growth_equation(
        const stat_t<real_t> &stat,
        const envi_t<real_t> &envi
      ) : stat(stat), envi(envi) {}

      // overloaded () operator invoked by thrust::transform()
      public: real_t operator()(const thrust_size_t id)
      {
        thrust_size_t ij = stat.ij[id]; 
        real_t drdt =
          ( // vapour density difference
            real_t(envi.rhod_rv[ij]) / si::cubic_metres * si::kilograms - // ambient rho_v
            ( // drop surface rho_v
              
                phc::p_vs<real_t>(envi.T[ij] * si::kelvins)
                / phc::R_v<real_t>() 
                / (envi.T[ij] * si::kelvins)
                * phc::a_w<real_t>(
                  this->rw3_of_xi(stat.xi[id]) * si::cubic_metres, // rw3
                  stat.rd3[id] * si::cubic_metres, // rd3
                  0 + stat.kpa[id] // kappa
                )
		// TODO: the Kelvin term
              
            )
          ) / phc::rho_w<real_t>() / (this->rw_of_xi(stat.xi[id]) * si::metres)
          * real_t(1e-5) * si::square_metres / si::seconds  // diffusion 
          / si::metres * si::seconds       // to make it dimensionless
          ; // TODO: move it to a phc_maxwell_mason.hpp
//cerr << "dxi_dt = " << this->dxidrw(this->rw_of_xi(stat.xi[id])) << " * " << drdt << "=" << this->dxidrw(this->rw_of_xi(stat.xi[id])) * drdt << endl;
        return this->dxidrw(this->rw_of_xi(stat.xi[id])) * drdt;// Jacobian
      }
    };

    // private fields
    private: stat_t<real_t> &stat;
    private: const envi_t<real_t> &envi;

    // ctor
    public: ode_xi(
      stat_t<real_t> &stat,
      const envi_t<real_t> &envi
    ) : stat(stat), envi(envi)
    {}

    // a post-ctor init method
    // (cannot be placed in the constructor as at that time the initial values were not loaded yet to psi)
    void init()
    {
      // initialising the wet radii (xi variable) to equilibrium values
      thrust::counting_iterator<thrust_size_t> iter(0);
      thrust::for_each(iter, iter + stat.n_part, equil(stat, envi));
    }
  
    // overloaded () operator invoked by odeint
    public: void operator()(
      const thrust::device_vector<real_t>&, 
      thrust::device_vector<real_t> &dxi_dt, 
      const real_t 
    )
    {
      thrust::counting_iterator<thrust_size_t> iter(0);
      thrust::transform(
        iter, iter + stat.n_part, 
        dxi_dt.begin(), 
        drop_growth_equation(stat, envi)
      );
    }
  };
}
