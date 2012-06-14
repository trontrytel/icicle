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
    static T dxidrw(const T &xi) { return 1; }
  };

  // xi = ln(rw / nm)
  template <typename T> struct xi_ln
  { 
    static T xi_of_rw(const T &rw) { return log(rw / T(1e-9)); }
    static T xi_of_rw3(const T &rw3) { return log(pow(rw3, T(1./3)) / T(1e-9)); }
    static T rw_of_xi(const T &xi) { return T(1e-9) * exp(xi); } 
    static T dxidrw(const T &rw) { return T(1) / rw; }
  };

  // xi = pow(rw, 2) 
  template <typename T> struct xi_p2
  { 
    static T xi_of_rw(const T &rw) { return rw * rw; }
    static T xi_of_rw3(const T &rw3) { return pow(rw3, T(2./3)); }
    static T rw_of_xi(const T &xi) { return sqrt(xi); } 
    static T dxidrw(const T &rw) { return T(2) * rw; }
  };

  // xi = pow(rw, 3) 
  template <typename T> struct xi_p3
  { 
    static T xi_of_rw(const T &rw) { return rw * rw * rw; }
    static T xi_of_rw3(const T &rw3) { return rw3; }
    static T rw_of_xi(const T &xi) { return pow(xi, T(1)/T(3)); } 
    static T dxidrw(const T &rw) { return T(3) * rw * rw; }
  };

  // nested the ODE system to be solved to update super-droplet sizes
  template <typename real_t, class algo, class xi>
  class ode_xi : public ode_algo<real_t, algo>
  { 
    public: real_t transform(const real_t &x) const
    {   
      return xi::xi_of_rw(x);
    }  

    // nested functor: setting wet radii to equilibrium values
    struct equil : xi
    { 
      stat_t<real_t> &stat;
      const envi_t<real_t> &envi;
      thrust::device_vector<real_t> &tmp;

      // ctor
      equil(
        stat_t<real_t> &stat, 
        const envi_t<real_t> &envi, 
        thrust::device_vector<real_t> &tmp
      ) : stat(stat), envi(envi), tmp(tmp) 
      {
        // nested functor
        struct init_tmp
        {
          const envi_t<real_t> &envi;
          thrust::device_vector<real_t> &tmp;
            
          // ctor
          init_tmp(const envi_t<real_t> &envi, thrust::device_vector<real_t> &tmp) 
            : envi(envi), tmp(tmp) 
          {}

          // overloded operator invoked by for_each below
          void operator()(thrust_size_t ij)
          {
            quantity<phc::mixing_ratio, real_t> 
              r = envi.rhod_rv[ij] / envi.rhod[ij];
            quantity<si::pressure, real_t> 
              p = phc::p<real_t>(envi.rhod_th[ij] * si::kelvins * si::kilograms / si::cubic_metres, r);
            quantity<si::temperature, real_t> 
              T = phc::T<real_t>((envi.rhod_th[ij] / envi.rhod[ij]) * si::kelvins, p, r);
  
            // - assuption of dr/dt = 0 => consistent only under subsaturation!
            // - initial values for vapour pressure field de facto assume aerosol at equilibrium
            // => using min(p/ps, .99) // TODO 99 as an option!
            tmp[ij] = thrust::min(
              real_t(.99),
              real_t(
                phc::R_v<real_t>() * T * (envi.rhod_rv[ij] * si::kilograms / si::cubic_metres) // TODO: interpolation to drop positions?
                / phc::p_vs<real_t>(T)
              )
            );
          } 
        };

        thrust::counting_iterator<thrust_size_t> iter(0);
        thrust::for_each(iter, iter + envi.n_cell, init_tmp(envi, tmp));
      }

      void operator()(thrust_size_t idx) 
      { 
        // ij and already sorted here!!!
        stat.xi[stat.id[idx]] = this->xi_of_rw3(phc::rw3_eq<real_t>(
          stat.rd3[stat.id[idx]] * si::cubic_metres, 
          0 + stat.kpa[stat.id[idx]], // it fails to compile without the zero!
          0 + tmp[stat.ij[idx]] // ditto
        ) / si::cubic_metres); 
      }
    };

    // nested functor: the drop growth law
    private: 
    class drop_growth_equation : public xi
    {
      // private fields
      private: const stat_t<real_t> &stat;

      // ctor
      public: drop_growth_equation(const stat_t<real_t> &stat) : stat(stat) {}

      // overloaded () operator invoked by thrust::transform()
      public: real_t operator()(thrust_size_t id)
      {
        return 0;
        //return this->dxidrw(this->rw(stat.xi[id]));
        //  * Maxwell-Mason;
      }
    };

    // private fields
    private: stat_t<real_t> &stat;
    private: const envi_t<real_t> &envi;
    thrust::device_vector<real_t> &tmp;

    // ctor
    public: ode_xi(
      stat_t<real_t> &stat,
      const envi_t<real_t> &envi,
      thrust::device_vector<real_t> &tmp
    ) : stat(stat), envi(envi), tmp(tmp)
    {}

    // a post-ctor init method
    // (cannot be placed in the constructor as at that time the initial values were not loaded yet to psi)
    void init()
    {
      // initialising the wet radii (xi variable) to equilibrium values
      thrust::counting_iterator<thrust_size_t> iter(0);
      thrust::for_each(iter, iter + stat.n_part, equil(stat, envi, tmp));
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
        drop_growth_equation(stat)
      );
    }
  };
}
