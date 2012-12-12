/** @file
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date July 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#pragma once
#include "sdm_ode.hpp"
#include "../phc/phc_henry.hpp"
#include "../phc/phc_chem_molar_mass.hpp"
#include "../phc/phc_dissociation.hpp" 
#include "../phc/phc_react.hpp"

namespace sdm
{
  template <typename real_t, class algo, class xi>
  class ode_chem : public ode_algo<real_t, algo>
  { 
    // nested functor 
    template <int chem, class nest_t>
    struct rhs : public xi
    {
      const nest_t &nest;

      // ctor
      rhs(const nest_t &nest) : nest(nest) {}

      // overloaded operator invoked by thrust
      real_t operator()(const thrust_size_t id)
      {
        switch (chem)
        {
          // Henrys law (see Seinfeld & Pandis 1997 p 340)
          // concentration of dissolved A = partial pressure of gas phase of A * Henrys coefficient for A
          case SO2: { 
            const thrust_size_t ij = nest.stat.ij[id];
            real_t c_equil = phc::henry::H_SO2<real_t>() * nest.opt_gas.at(gSO2) * (nest.envi.p[ij] * si::pascals)  
              * (4./3 * M_PI * this->rw3_of_xi(nest.stat.xi[id]) * si::cubic_metres) * phc::mass::M_SO2<real_t>() / si::kilograms;
            return (c_equil - nest.stat.c_aq[id + nest.stat.n_part  * SO2]) / nest.dt * si::seconds;
          }
          case O3: {
            const thrust_size_t ij = nest.stat.ij[id];
            real_t c_equil = phc::henry::H_O3<real_t>() * nest.opt_gas.at(gO3) * (nest.envi.p[ij] * si::pascals)  
              * (4./3 * M_PI * this->rw3_of_xi(nest.stat.xi[id]) * si::cubic_metres) * phc::mass::M_O3<real_t>() / si::kilograms;
            return (c_equil - nest.stat.c_aq[id + nest.stat.n_part * O3]) / nest.dt * si::seconds;
          }
          case H2O2: {
            const thrust_size_t ij = nest.stat.ij[id];
            real_t c_equil = phc::henry::H_H2O2<real_t>() * nest.opt_gas.at(gH2O2) * (nest.envi.p[ij] * si::pascals)  
              * (4./3 * M_PI * this->rw3_of_xi(nest.stat.xi[id]) * si::cubic_metres) * phc::mass::M_H2O2<real_t>() / si::kilograms;
            return (c_equil - nest.stat.c_aq[id + nest.stat.n_part * H2O2]) / nest.dt * si::seconds;
          }

          // dissociation (see Seinfeld & Pandis 1997 p 349)
          // [A-] = dissociation const * [H2O * A] / [H+]
          case H: {
/*          const thrust_size_t ij = nest.stat.ij[id];
            real_t c_equil = (
               nest.stat.c_aq[id + nest.stat.n_part * OH]   / (phc::mass::M_OH<real_t>() * si::moles / si::kilograms)
            +  nest.stat.c_aq[id + nest.stat.n_part * HSO3] / (phc::mass::M_HSO3<real_t>() * si::moles / si::kilograms)
//             + nest.stat.c_aq[id + nest.stat.n_part * SO3]  / (phc::mass::M_SO3<real_t>() * si::moles / si::kilograms)
            )* phc::mass::M_H<real_t>() * si::moles / si::kilograms ;
            return (c_equil - nest.stat.c_aq[id + nest.stat.n_part * H]) / nest.dt * si::seconds;
*/
return 0;
          }
          case OH: {
            const thrust_size_t ij = nest.stat.ij[id];                         
            real_t c_equil = 
             (4./3 * M_PI * this->rw3_of_xi(nest.stat.xi[id]) * si::cubic_metres)
             * (4./3 * M_PI * this->rw3_of_xi(nest.stat.xi[id]) * si::cubic_metres)
             * phc::dissociation::K_H2O<real_t>()
             * phc::mass::M_H<real_t>() / si::kilograms
             * phc::mass::M_OH<real_t>() / si::kilograms
             / real_t(nest.stat.c_aq[id + nest.stat.n_part * H]);
           return (c_equil - nest.stat.c_aq[id + nest.stat.n_part * OH]) / nest.dt * si::seconds;
          }

          case HSO3: {
            const thrust_size_t ij = nest.stat.ij[id];                         
            real_t c_equil = (4./3 * M_PI * this->rw3_of_xi(nest.stat.xi[id]) * si::cubic_metres) 
              * phc::mass::M_HSO3<real_t>() * phc::dissociation::K_SO2<real_t>() / si::kilograms
              * real_t(nest.stat.c_aq[id + nest.stat.n_part * SO2]) / real_t(nest.stat.c_aq[id + nest.stat.n_part * H]) 
              * (phc::mass::M_H<real_t>() / phc::mass::M_SO2<real_t>());
            return (c_equil - nest.stat.c_aq[id + nest.stat.n_part * HSO3]) / nest.dt * si::seconds;
          }
          case SO3: {
            const thrust_size_t ij = nest.stat.ij[id];                         
/*
            real_t c_equil = (4./3 * M_PI * this->rw3_of_xi(nest.stat.xi[id]) * si::cubic_metres) 
              * phc::mass::M_SO3<real_t>() * phc::dissociation::K_HSO3<real_t>() / si::kilograms 
              * real_t(nest.stat.c_aq[id + nest.stat.n_part * HSO3]) / real_t(nest.stat.c_aq[id + nest.stat.n_part * H]) 
              * (phc::mass::M_H<real_t>() / phc::mass::M_HSO3<real_t>()) ;
            return (c_equil - nest.stat.c_aq[id + nest.stat.n_part * SO3]) / nest.dt * si::seconds;
*/
return 0;
          }

          case HSO4: {
            const thrust_size_t ij = nest.stat.ij[id];                   
/*                       // should be si::mass_density ? 
            quantity<si::concentration, real_t> dS_VI = 
              real_t(nest.stat.c_aq[id + nest.stat.n_part * O3]) * si::kilograms / phc::mass::M_O3<real_t>() 
              / (4./3 * M_PI * this->rw3_of_xi(nest.stat.xi[id]) * si::cubic_metres)  
              * (
                real_t(nest.stat.c_aq[id + nest.stat.n_part * SO2]) / phc::mass::M_SO2<real_t>() * phc::react::R_S_O3_k0<real_t>() +  
                real_t(nest.stat.c_aq[id + nest.stat.n_part * HSO3]) / phc::mass::M_HSO3<real_t>() * phc::react::R_S_O3_k1<real_t>() +
                real_t(nest.stat.c_aq[id + nest.stat.n_part * SO3]) / phc::mass::M_SO3<real_t>() * phc::react::R_S_O3_k2<real_t>() 
              ) * si::kilograms / (4./3 * M_PI * this->rw3_of_xi(nest.stat.xi[id]) * si::cubic_metres) 
              * nest.dt;

            quantity<si::concentration, real_t> H = real_t(nest.stat.c_aq[id + nest.stat.n_part * H]) * si::kilograms 
              / phc::mass::M_H<real_t>() / (4./3 * M_PI * this->rw3_of_xi(nest.stat.xi[id]) * si::cubic_metres)
            return H * (dS_VI / (H + phc::dissociation::K_HSO4<real_t>());
*/
return 0;
          }
          case SO4: {
            return 0;
          }
        }
      }
    };

    // private fields
    private: const stat_t<real_t> &stat;
    private: const envi_t<real_t> &envi;
    private: const thrust::counting_iterator<thrust_size_t> zero;
    private: quantity<si::time, real_t> dt;
    private: const map<enum chem_gas, quantity<phc::mixing_ratio, real_t>> opt_gas;
    private: const map<enum chem_aq, quantity<si::mass, real_t>> opt_aq;

    void init(const quantity<si::time, real_t> &dt_arg)
    {
      dt = dt_arg;
    }

    // ctor 
    public: ode_chem(
      const stat_t<real_t> &stat, 
      const envi_t<real_t> &envi, 
      const map<enum chem_gas, quantity<phc::mixing_ratio, real_t>> &opt_gas,
      const map<enum chem_aq, quantity<si::mass, real_t>> opt_aq
    ) 
      : stat(stat), envi(envi), zero(0), opt_gas(opt_gas), opt_aq(opt_aq)
    {}

    // overloaded () operator invoked by odeint
    public: void operator()(
      const thrust::device_vector<real_t> &, 
      thrust::device_vector<real_t> &dc_dt, 
      const real_t
    )
    {
      thrust::transform(zero, zero + stat.n_part, dc_dt.begin() + SO2  * stat.n_part, 
        rhs<SO2, decltype(*this)>(*this)
      );
      thrust::transform(zero, zero + stat.n_part, dc_dt.begin() + O3 * stat.n_part, 
        rhs<O3, decltype(*this)>(*this)
      );
      thrust::transform(zero, zero + stat.n_part, dc_dt.begin() + H2O2 * stat.n_part, 
        rhs<H2O2, decltype(*this)>(*this)
      );
      thrust::transform(zero, zero + stat.n_part, dc_dt.begin() + H * stat.n_part, 
        rhs<H, decltype(*this)>(*this)
      );
      thrust::transform(zero, zero + stat.n_part, dc_dt.begin() + OH * stat.n_part, 
        rhs<OH, decltype(*this)>(*this)
      );
      thrust::transform(zero, zero + stat.n_part, dc_dt.begin() + HSO3 * stat.n_part, 
        rhs<HSO3, decltype(*this)>(*this)
      );
      thrust::transform(zero, zero + stat.n_part, dc_dt.begin() + SO3 * stat.n_part, 
        rhs<SO3, decltype(*this)>(*this)
      );
      thrust::transform(zero, zero + stat.n_part, dc_dt.begin() + HSO4 * stat.n_part, 
        rhs<HSO4, decltype(*this)>(*this)
      );
      thrust::transform(zero, zero + stat.n_part, dc_dt.begin() + SO4 * stat.n_part, 
        rhs<SO4, decltype(*this)>(*this)
      );
    }
  };
}
