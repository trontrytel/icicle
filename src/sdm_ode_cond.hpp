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
  class ode_cond : public ode_algo<real_t, algo>
  { 
    // nested functor
    class m_3 : xi
    {
      private: real_t dv;
      private: const stat_t<real_t> &stat;
      private: thrust::device_vector<real_t> &out;
 
      public: m_3(
        const stat_t<real_t> &stat,
        thrust::device_vector<real_t> &out,
        const grd<real_t> &grid
      ) :
        stat(stat),
        out(out),
        dv(grid.dx() * grid.dy() * grid.dz() / si::cubic_metres)
      {
        thrust::fill(out.begin(), out.end(), real_t(0));
      }
    
      void operator()(const thrust_size_t id)
      {
        out[stat.ij[id]] += real_t(stat.n[id]) * this->rw3_of_xi(stat.xi[id]) / dv;
//cerr << "!!!" << stat.n[id] << " " << stat.xi[id] << " " << this->rw3_of_xi(stat.xi[id]) << endl;
      }
    };

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

      void operator()(const thrust_size_t id) 
      { 
        thrust_size_t ij = stat.ij[id]; 
        real_t rw3_eq = phc::rw3_eq<real_t>( // TODO: allow choice among different Koehler curves
          stat.rd3[id] * si::cubic_metres, 
          0 + stat.kpa[id], // it fails to compile without the zero! // TODO: real_t()
          thrust::min(
            real_t(.99),
            real_t( // TODO: interpolation to drop positions?
              phc::R_v<real_t>() * (envi.T[ij] * si::kelvins) * (envi.rhod_rv[ij] * si::kilograms / si::cubic_metres) 
              / (phc::p_vs<real_t>(envi.T[ij] * si::kelvins))
            )
          )
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

      // ctor
      public: drop_growth_equation(
        const stat_t<real_t> &stat,
        const envi_t<real_t> &envi,
        const real_t dt
      ) : stat(stat), envi(envi), dt(dt) {}

      // overloaded () operator invoked by thrust::transform()
      public: real_t operator()(const thrust_size_t id)
      {
        thrust_size_t ij = stat.ij[id]; 
        real_t rdrdt =
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
          ) / phc::rho_w<real_t>() 
          * real_t(2.21e-5) * si::square_metres / si::seconds  // diffusion 
          * si::seconds / si::square_metres       // to make it dimensionless
          ; // TODO: move it to a phc_maxwell_mason.hpp
        real_t dxidt = rdrdt * this->dxidrw(this->rw_of_xi(stat.xi[id])) / this->rw_of_xi(stat.xi[id]);
return std::max(real_t(0), dxidt); // TODO: only if evap turned off
/*
        if (dt != 0 && stat.xi[id] + dt * dxidt < 0) dxidt = (this->xi_of_rw3(
          phc::rw3_eq<real_t>(
            stat.rd3[id] * si::cubic_metres, 
            real_t(stat.kpa[id]),
            real_t( // TODO: interpolation to drop positions?
              thrust::min(
                real_t(.99),
                real_t(phc::R_v<real_t>() * (envi.T[ij] * si::kelvins) * (envi.rhod_rv[ij] * si::kilograms / si::cubic_metres) 
                / (phc::p_vs<real_t>(envi.T[ij] * si::kelvins)))
              )   
            )
          ) / si::cubic_metres) 
          - stat.xi[id]
        ) / dt;
if (stat.xi[id] < 0 || !isfinite(dxidt)) 
{
  cerr << "dxidt = " << rdrdt << " * " << this->dxidrw(this->rw_of_xi(stat.xi[id])) / this->rw_of_xi(stat.xi[id]) << endl;
  cerr << "xi_new = xi + dt * dxidt = " << stat.xi[id] << " + " << dt << " * " << dxidt << endl;
  exit(1);
}
        return .1 *  dxidt;
*/
// TODO? this->dxidrw_by_r?
      }
    };

    // private fields
    private: stat_t<real_t> &stat;
    private: envi_t<real_t> &envi;
    private: const grd<real_t> &grid;

    // ctor
    public: ode_cond(
      stat_t<real_t> &stat,
      envi_t<real_t> &envi,
      const grd<real_t> &grid
    ) : stat(stat), envi(envi), grid(grid)
    {}
  
    // overloaded () operator invoked by odeint
    public: void operator()(
      const thrust::device_vector<real_t>&, 
      thrust::device_vector<real_t> &dxi_dt, 
      const real_t t
    )
    {
      thrust::counting_iterator<thrust_size_t> iter(0);
      thrust::transform(
        iter, iter + stat.n_part, 
        dxi_dt.begin(), 
        drop_growth_equation(stat, envi, t) // here: t = dt !
      );
    }

    // a post-ctor init method
    // (cannot be placed in the constructor as at that time the initial values were not loaded yet to psi)
    void init()
    {
      thrust::counting_iterator<thrust_size_t> iter(0);

      // initialising the wet radii (xi variable) to equilibrium values
      thrust::for_each(iter, iter + stat.n_part, equil(stat, envi));

      // 
      thrust::for_each(iter, iter + stat.n_part, m_3(stat, envi.m_3_old, grid));
    }

    // a post-ode-step method
    void adjust()
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
            //if (rhod_rv != rhod_rv_initial) // TODO 
            //{
              envi.r[ij] = rhod_rv / (envi.rhod[ij] * si::kilograms / si::cubic_metres);
              envi.p[ij] = phc::p<real_t>(rhod_th, real_t(envi.r[ij])) / si::pascals;
              envi.T[ij] = phc::T<real_t>(
                rhod_th / (envi.rhod[ij] * si::kilograms / si::cubic_metres), 
                envi.p[ij] * si::pascals, 
                real_t(envi.r[ij])
              ) / si::kelvins;
            //}
            F = phc::dtheta_drv<real_t>(
              envi.T[ij] * si::kelvins, 
              envi.p[ij] * si::pascals, 
              real_t(envi.r[ij]), 
              rhod_th, 
              envi.rhod[ij] * si::kilograms / si::cubic_metres
            ); // TODO: option : which dtheta...
          }
        };

        // odeint::euler< // TODO: opcja?
        odeint::runge_kutta4<
          quantity<multiply_typeof_helper<si::mass_density, si::temperature>::type, real_t>, // state_type
          real_t, // value_type
          quantity<si::temperature, real_t>, // deriv_type
          quantity<si::mass_density, real_t>, // time_type
          odeint::vector_space_algebra, 
          odeint::default_operations, 
          odeint::never_resizer
        > S; // TODO: instantiate it in the ctor?

        private: envi_t<real_t> &envi;

        // ctor
        public: adj(envi_t<real_t> &envi) : envi(envi) { }

        // invoked by thrust::transform
        public: void operator()(const thrust_size_t ij) // TODO: too much duplication with eqs_todo_bulk_ode.cpp :(
        {
          // theta is modified by do_step, and hence we cannot pass an expression and we need a temp. var.
          quantity<multiply_typeof_helper<si::mass_density, si::temperature>::type, real_t>
            tmp = envi.rhod_th[ij] * si::kilograms / si::cubic_metres * si::kelvins;

          // (<r^3>_new - <r^3>_old)  --->  drhod_rv
          quantity<si::mass_density, real_t> drhod_rv =
            (envi.m_3_old[ij] - envi.m_3_new[ij]) * phc::rho_w<real_t>() * real_t(4./3) * phc::pi<real_t>();

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

          // updating rhod_rv
          envi.rhod_rv[ij] += drhod_rv / (si::kilograms / si::cubic_metres);
          //assert(rhod_rv(i,j,k) >= 0);  TODO
          //assert(isfinite(rhod_rv(i,j,k))); TODO
        }
      };
   
      thrust::counting_iterator<thrust_size_t> zero(0);

      // calculating new <r^3>
      thrust::for_each(zero, zero + stat.n_part, m_3(stat, envi.m_3_new, grid)); // TODO: does the zero work????

      // adjusting rhod_rv and rhod_th according to the First Law (and T,p,r as well)
      thrust::for_each(zero, zero + envi.n_cell, adj(envi));

      // new -> old 
      thrust::swap(envi.m_3_old, envi.m_3_new); 
    }
  };
} // namespace sdm
