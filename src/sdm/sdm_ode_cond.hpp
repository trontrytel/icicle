/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date March-June 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#pragma once

#include "sdm_ode.hpp"
#include "../phc/phc_theta.hpp"
#include "../phc/phc_const_cp.hpp"
#include "../phc/phc_kappa_koehler.hpp"
#include "../phc/phc_kelvin_term.hpp"
#include "../phc/phc_maxwell-mason.hpp"

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
  class ode_cond : public ode_algo<real_t, algo>
  { 
    // nested functor
    class dm_3_summator : xi
    {
      private: const real_t dv, sign;
      private: const stat_t<real_t> &stat;
      private: thrust::device_vector<real_t> &out;
 
      // ctor
      public: dm_3_summator(
        const stat_t<real_t> &stat,
        thrust::device_vector<real_t> &out,
        const grd<real_t> &grid,
        const real_t sign
      ) :
        stat(stat),
        out(out),
        dv(grid.dx() * grid.dy() * grid.dz() / si::cubic_metres),
        sign(sign)
      { }
    
      void operator()(const thrust_size_t id)
      {
        out[stat.ij[id]] += sign * real_t(stat.n[id]) * this->rw3_of_xi(stat.xi[id]) / dv;
      }
    };

    // nested functor: setting wet radii to equilibrium values
    struct equil : xi
    { 
      stat_t<real_t> &stat;
      const envi_t<real_t> &envi;
      const real_t maxRH;

      // ctor
      equil(
        stat_t<real_t> &stat, 
        const envi_t<real_t> &envi,
        const real_t maxRH
      ) : stat(stat), envi(envi), maxRH(maxRH)
      { }

      void operator()(const thrust_size_t id) 
      { 
        thrust_size_t ij = stat.ij[id]; 
        real_t rw3_eq = phc::kappa::rw3_eq<real_t>( // TODO: allow choice among different Koehler curves
          stat.rd3[id] * si::cubic_metres, 
          real_t(stat.kpa[id]),
          thrust::min(
            maxRH,
            real_t( // TODO: interpolation to drop positions?
              phc::R_v<real_t>() * (envi.T[ij] * si::kelvins) * (envi.rhod_rv[ij] * si::kilograms / si::cubic_metres) 
              / (phc::p_vs<real_t>(envi.T[ij] * si::kelvins))
            )
          ),
          envi.T[ij] * si::kelvins
        ) / si::cubic_metres;
        stat.xi[id] = this->xi_of_rw3(rw3_eq);
      }
    };

    // nested functor: the drop growth law
    private: 
    class drop_growth_equation : public xi
    {
      // private fields
      private: const stat_t<real_t> &stat;
      private: const envi_t<real_t> &envi;
      private: const real_t dt;
      private: const real_t maxRH;

      // ctor
      public: drop_growth_equation(
        const stat_t<real_t> &stat,
        const envi_t<real_t> &envi,
        const real_t dt,
        const real_t maxRH
      ) : stat(stat), envi(envi), dt(dt), maxRH(maxRH) {}

      // overloaded () operator invoked by thrust::transform()
      public: real_t operator()(const thrust_size_t id)
      {
        thrust_size_t ij = stat.ij[id]; 

        // the growth rate defined by the Maxwell-Mason equation
        real_t dxidt = phc::rdrdt<real_t>(
          real_t(envi.rhod_rv[ij]) / si::cubic_metres * si::kilograms,
          envi.T[ij] * si::kelvins,
          phc::kappa::a_w<real_t>(
            this->rw3_of_xi(stat.xi[id]) * si::cubic_metres, // rw3
            stat.rd3[id] * si::cubic_metres, // rd3
            real_t(stat.kpa[id]) // kappa
          ),
          phc::kelvin::klvntrm<real_t>(this->rw_of_xi(stat.xi[id]) * si::metres, envi.T[ij] * si::kelvins) // the Kelvin term
        ) 
        * si::seconds / si::square_metres // to make it dimensionless
        / this->rw_of_xi(stat.xi[id]) // to get drdt
        * this->dxidrw(this->rw_of_xi(stat.xi[id])) // to get dxidt (Jacobian)
        ;

        // the growth rate defined by "come back to equilibrium"
        // compare e.g. apparently analogous tricks in:
        // - Johnson 1980, JAS 37 2079-2085, disussion of inequality (5)
        // - Heanel 1987, Beitr. Phys. Atmosph. 60 321-338, section 4.2
        real_t xi_eq = this->xi_of_rw3(
          phc::kappa::rw3_eq<real_t>(
            stat.rd3[id] * si::cubic_metres, 
            real_t(stat.kpa[id]),
            real_t( // TODO: interpolation to drop positions?
              thrust::min(
                maxRH,
                real_t(phc::R_v<real_t>() * (envi.T[ij] * si::kelvins) * (envi.rhod_rv[ij] * si::kilograms / si::cubic_metres) 
                / (phc::p_vs<real_t>(envi.T[ij] * si::kelvins)))
              )   
            ),
            envi.T[ij] * si::kelvins
          ) / si::cubic_metres);

        // TODO: as an option
        // chosing which one to use
        if (stat.xi[id] + dt * dxidt < xi_eq) dxidt = (xi_eq - stat.xi[id]) / dt;

        // returning
        return dxidt;
// TODO? this->dxidrw_by_r?
      }
    };

    // private fields
    private: stat_t<real_t> &stat;
    private: envi_t<real_t> &envi;
    private: const grd<real_t> &grid;
    private: quantity<si::time, real_t> dt_sstp;
    private: const real_t maxRH = .99; // TODO: should this be an option (at least it deserves some documentation)
    private: thrust::device_vector<real_t> &dm_3;

    // ctor
    public: ode_cond(
      stat_t<real_t> &stat,
      envi_t<real_t> &envi,
      const grd<real_t> &grid,
      thrust::device_vector<real_t> &tmp_real
    ) : stat(stat), envi(envi), grid(grid), dm_3(tmp_real)
    {}
  
    // a post-ctor init method
    // (cannot be placed in the constructor as at that time the initial values were not loaded yet to psi)
    void init(const quantity<si::time, real_t> &dt_sstp_)
    {
      dt_sstp = dt_sstp_;      

      // initialising the wet radii (xi variable) to equilibrium values
      thrust::counting_iterator<thrust_size_t> zero(0);
      thrust::for_each(zero, zero + stat.n_part, equil(stat, envi, maxRH));
    }

    //
    void pre_step()
    {
      thrust::counting_iterator<thrust_size_t> zero(0);

      // calculating the third moment of the spectrum before condensation
      thrust::fill(dm_3.begin(), dm_3.end(), real_t(0));
      thrust::for_each(zero, zero + stat.n_part, dm_3_summator(stat, dm_3, grid, -1));
    }

    // overloaded () operator invoked by odeint
    public: void operator()(
      const thrust::device_vector<real_t>&, 
      thrust::device_vector<real_t> &dxi_dt, 
      const real_t 
    )
    {
      thrust::counting_iterator<thrust_size_t> zero(0);

      // calculating drop growth
      thrust::transform(
        zero, zero + stat.n_part, 
        dxi_dt.begin(), 
        drop_growth_equation(stat, envi, dt_sstp / si::seconds, maxRH)
      );
    }

    //
    void post_step()
    {
      // nested functor
      class adj
      {
        // nested functor 
        class rhs
        {
          private: const thrust_size_t ij;
          private: envi_t<real_t> &envi;
 
          // ctor
          public: rhs(
            envi_t<real_t> &envi,
            const thrust_size_t ij
          ) : envi(envi), ij(ij) {}
           
          // overloaded op. invoked by odeint
          public: void operator()(
            const quantity<multiply_typeof_helper<si::mass_density, si::temperature>::type, real_t> rhod_th,
            quantity<si::temperature, real_t> &F,
            const quantity<si::mass_density, real_t> rhod_rv
          )
          {
            // TODO: if other than Euler, the envi.{p,T,r} should be updated here!
            F = phc::dtheta_drv<real_t>(
              envi.T[ij] * si::kelvins, 
              envi.p[ij] * si::pascals, 
              real_t(envi.r[ij]), 
              rhod_th, 
              envi.rhod[ij] * si::kilograms / si::cubic_metres
            ); 
          }
        };

        // odeint::euler< // TODO: opcja?
        odeint::euler<
          quantity<multiply_typeof_helper<si::mass_density, si::temperature>::type, real_t>, // state_type
          real_t, // value_type
          quantity<si::temperature, real_t>, // deriv_type
          quantity<si::mass_density, real_t>, // time_type
          odeint::vector_space_algebra, 
          odeint::default_operations, 
          odeint::never_resizer
        > S; // TODO: instantiate it in the ctor?

        private: envi_t<real_t> &envi;
        private: const thrust::device_vector<real_t> &dm_3;

        // ctor
        public: adj(
          envi_t<real_t> &envi,
          thrust::device_vector<real_t> &dm_3
        ) : envi(envi), dm_3(dm_3) 
        { }

        // invoked by thrust::transform
        public: void operator()(const thrust_size_t ij) // TODO: too much duplication with eqs_todo_bulk_ode.cpp :(
        {
          // (<r^3>_new - <r^3>_old)  --->  drhod_rv
          quantity<si::mass_density, real_t> drhod_rv =  
            real_t(-dm_3[ij]) * phc::rho_w<real_t>() * real_t(4./3) * phc::pi<real_t>();

          {
            // theta is modified by do_step, and hence we cannot pass an expression and we need a temp. var.
            quantity<multiply_typeof_helper<si::mass_density, si::temperature>::type, real_t>
              tmp = envi.rhod_th[ij] * si::kilograms / si::cubic_metres * si::kelvins;

            // integrating the First Law for moist air
            rhs F(envi, ij); // TODO: instantiate it somewehere else - not every sub-time-step !!!
            S.do_step(
              boost::ref(F),
              tmp,
              envi.rhod_rv[ij] * si::kilograms / si::cubic_metres,
              drhod_rv 
            );

            // latent heat source/sink due to evaporation/condensation
            envi.rhod_th[ij] = tmp / (si::kilograms / si::cubic_metres * si::kelvins); 
          }

          // updating rhod_rv
          envi.rhod_rv[ij] += drhod_rv / (si::kilograms / si::cubic_metres);
          {
            envi.r[ij] = envi.rhod_rv[ij] / envi.rhod[ij];
            envi.p[ij] = phc::p<real_t>(envi.rhod_th[ij] * si::kelvins * si::kilograms / si::cubic_metres, real_t(envi.r[ij])) / si::pascals;
            envi.T[ij] = phc::T<real_t>(
              envi.rhod_th[ij] / envi.rhod[ij] * si::kelvins, 
              envi.p[ij] * si::pascals, 
              real_t(envi.r[ij])
            ) / si::kelvins;
          }
          //assert(rhod_rv(i,j,k) >= 0);  TODO
          //assert(isfinite(rhod_rv(i,j,k))); TODO
        }
      };
   
      thrust::counting_iterator<thrust_size_t> zero(0);

      // substracting the third moment after condensation
      thrust::for_each(zero, zero + stat.n_part, dm_3_summator(stat, dm_3, grid, +1));

      // adjusting rhod_rv and rhod_th according to the First Law (and T,p,r as well)
      thrust::for_each(zero, zero + envi.n_cell, adj(envi, dm_3));
    }
  };
} // namespace sdm
