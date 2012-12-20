/** @file
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date July 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @brief definition of source terms for aqueous phase reactions of oxidation of S(IV) by O3 and H2O2 
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
      const thrust_size_t ij = nest.stat.ij[id];                   
      quantity<si::volume, real_t> V = 4./3 * M_PI * this->rw3_of_xi(nest.stat.xi[id]) * si::cubic_metres; //volume of the droplet

      // helpers for O3 reactions
      quantity<si::mass, real_t> O3_SO2 = real_t(nest.stat.c_aq[id + nest.stat.n_part * O3]) * si::kilograms / V * nest.dt
        * real_t(nest.stat.c_aq[id + nest.stat.n_part * SO2]) * si::kilograms / phc::mass::M_SO2<real_t>() * phc::react::R_S_O3_k0<real_t>();  
      quantity<si::mass, real_t> O3_HSO3 = real_t(nest.stat.c_aq[id + nest.stat.n_part * O3]) * si::kilograms / V * nest.dt
        * real_t(nest.stat.c_aq[id + nest.stat.n_part * HSO3]) * si::kilograms / phc::mass::M_HSO3<real_t>() * phc::react::R_S_O3_k1<real_t>();
      quantity<si::mass, real_t> O3_SO3 = real_t(nest.stat.c_aq[id + nest.stat.n_part * O3]) * si::kilograms / V * nest.dt
        * real_t(nest.stat.c_aq[id + nest.stat.n_part * SO3]) * si::kilograms / phc::mass::M_SO3<real_t>() * phc::react::R_S_O3_k2<real_t>(); 

      // helper for H2O2 reactions
      quantity<si::amount, real_t> H2O2_HSO3 = phc::react::R_S_H2O2_k<real_t>() / pow<2>(V)
        / phc::mass::M_H2O2<real_t>() / phc::mass::M_H<real_t>() / phc::mass::M_HSO3<real_t>()
        * real_t(nest.stat.c_aq[id + nest.stat.n_part * H2O2]) * si::kilograms 
        * real_t(nest.stat.c_aq[id + nest.stat.n_part * H])    * si::kilograms 
        * real_t(nest.stat.c_aq[id + nest.stat.n_part * HSO3]) * si::kilograms / ( 
        1. + phc::react::R_S_H2O2_K<real_t>() 
             * real_t(nest.stat.c_aq[id + nest.stat.n_part * H]) * si::kilograms / phc::mass::M_H<real_t>() / V
        ) * nest.dt;

        switch (chem)
        {
          //S(IV)->S(VI) oxidation by H2O2 and O3
          case S_VI: {
            quantity<si::mass, real_t> dS_VI_O3 = phc::mass::M_H2SO4<real_t>() / phc::mass::M_O3<real_t>() * (O3_SO2 + O3_HSO3 + O3_SO3);
            quantity<si::mass, real_t> dS_VI_H2O2 = phc::mass::M_H2SO4<real_t>() * H2O2_HSO3; 

            return (dS_VI_O3 + dS_VI_H2O2) / si::kilograms;
          }
          case H2O2 : return - phc::mass::M_H2O2<real_t>() * H2O2_HSO3 / si::kilograms; 
          case O3 : return -(O3_SO2 + O3_HSO3 + O3_SO3) / si::kilograms;
          case SO2 : return - phc::mass::M_SO2<real_t>() / phc::mass::M_O3<real_t>() * O3_SO2 / si::kilograms;  
          case HSO3:{
            quantity<si::mass, real_t> dHSO3_O3 = phc::mass::M_HSO3<real_t>() / phc::mass::M_O3<real_t>() * O3_HSO3;
            quantity<si::mass, real_t> dHSO3_H2O2 = phc::mass::M_HSO3<real_t>() * H2O2_HSO3;
            
            return - (dHSO3_O3 + dHSO3_H2O2) / si::kilograms;
          }
          case SO3 : return - phc::mass::M_SO3<real_t>() / phc::mass::M_O3<real_t>() * O3_SO3 / si::kilograms; 
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
      thrust::transform(zero, zero + stat.n_part, dc_dt.begin() + S_VI * stat.n_part, rhs<S_VI, decltype(*this)>(*this));
      thrust::transform(zero, zero + stat.n_part, dc_dt.begin() + O3   * stat.n_part, rhs<O3,   decltype(*this)>(*this));
      thrust::transform(zero, zero + stat.n_part, dc_dt.begin() + H2O2 * stat.n_part, rhs<H2O2, decltype(*this)>(*this));
      thrust::transform(zero, zero + stat.n_part, dc_dt.begin() + SO2  * stat.n_part, rhs<SO2,  decltype(*this)>(*this));
      thrust::transform(zero, zero + stat.n_part, dc_dt.begin() + HSO3 * stat.n_part, rhs<HSO3, decltype(*this)>(*this));
      thrust::transform(zero, zero + stat.n_part, dc_dt.begin() + SO3  * stat.n_part, rhs<SO3,  decltype(*this)>(*this));
    }
  };
}
