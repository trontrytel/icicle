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
#include "phc_henry.hpp"
#include "phc_chem_molar_mass.hpp"
#include "phc_dissociation.hpp" 

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
          case sdm::SO2: { 
            const thrust_size_t ij = nest.stat.ij[id];
            real_t c_equil = phc::henry::H_SO2<real_t>() * nest.opt_gas.at(sdm::gSO2) * (nest.envi.p[ij] * si::pascals)  
              * si::cubic_metres / si::moles;
            return (c_equil - nest.stat.c_aq[id + nest.stat.n_part  * sdm::SO2]) / nest.dt * si::seconds;
          }
          case sdm::O3: {
            const thrust_size_t ij = nest.stat.ij[id];
            real_t c_equil = phc::henry::H_O3<real_t>() * nest.opt_gas.at(sdm::gO3) * (nest.envi.p[ij] * si::pascals)  
              * si::cubic_metres / si::moles;
            return (c_equil - nest.stat.c_aq[id + nest.stat.n_part * sdm::O3]) / nest.dt * si::seconds;
          }
          case sdm::H2O2: {
            const thrust_size_t ij = nest.stat.ij[id];
            real_t c_equil = phc::henry::H_H2O2<real_t>() * nest.opt_gas.at(sdm::gH2O2) * (nest.envi.p[ij] * si::pascals)  
              * si::cubic_metres / si::moles;
            return (c_equil - nest.stat.c_aq[id + nest.stat.n_part * sdm::H2O2]) / nest.dt * si::seconds;
          }
          // dissociation (see Seinfeld & Pandis 1997 p 349)
          // [A-] = dissociation const * [H2O * A] / [H+]
          case sdm::H: {
          const thrust_size_t ij = nest.stat.ij[id];
            real_t c_equil = 
               nest.stat.c_aq[id + nest.stat.n_part * sdm::OH] 
             + nest.stat.c_aq[id + nest.stat.n_part * sdm::HSO3]
             + nest.stat.c_aq[id + nest.stat.n_part * sdm::SO3];
            return (c_equil - nest.stat.c_aq[id + nest.stat.n_part * sdm::HSO3]) / nest.dt * si::seconds;
          }
          case sdm::OH: {
            const thrust_size_t ij = nest.stat.ij[id];                         
            real_t c_equil = phc::dissociation::K_H2O<real_t>() * si::cubic_metres / si::moles
             * real_t(1/18.*1e3) // [H2O]   TODO !!!!
             / real_t(nest.stat.c_aq[id + nest.stat.n_part * sdm::H]);
            return (c_equil - nest.stat.c_aq[id + nest.stat.n_part * sdm::OH]) / nest.dt * si::seconds;
          }
          case sdm::HSO3: {
            const thrust_size_t ij = nest.stat.ij[id];                         
            real_t c_equil = phc::dissociation::K_SO2<real_t>() * si::cubic_metres / si::moles
             * real_t(nest.stat.c_aq[id + nest.stat.n_part * sdm::SO2])
             / real_t(nest.stat.c_aq[id + nest.stat.n_part * sdm::H]);
            return (c_equil - nest.stat.c_aq[id + nest.stat.n_part * sdm::HSO3]) / nest.dt * si::seconds;
          }
          case sdm::SO3: {
            const thrust_size_t ij = nest.stat.ij[id];                         
            real_t c_equil = phc::dissociation::K_HSO3<real_t>() * si::cubic_metres / si::moles
             * real_t(nest.stat.c_aq[id + nest.stat.n_part * sdm::HSO3])
             / real_t(nest.stat.c_aq[id + nest.stat.n_part * sdm::H]);
            return (c_equil - nest.stat.c_aq[id + nest.stat.n_part * sdm::SO3]) / nest.dt * si::seconds;
          }
        }
      }
    };

    // private fields
    private: const stat_t<real_t> &stat;
    private: const envi_t<real_t> &envi;
    private: const thrust::counting_iterator<thrust_size_t> zero;
    private: quantity<si::time, real_t> dt;
    private: const map<enum sdm::chem_gas, quantity<phc::mixing_ratio, real_t>> opt_gas;
    private: const map<enum sdm::chem_aq, quantity<divide_typeof_helper<si::amount, si::volume>::type, real_t>> opt_aq;

    void init(const quantity<si::time, real_t> &dt_arg)
    {
      dt = dt_arg;
    }

    // ctor 
    public: ode_chem(
      const stat_t<real_t> &stat, 
      const envi_t<real_t> &envi, 
      const map<enum sdm::chem_gas, quantity<phc::mixing_ratio, real_t>> &opt_gas,
      const map<enum sdm::chem_aq, quantity<divide_typeof_helper<si::amount, si::volume>::type, real_t>> opt_aq
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
      thrust::transform(zero, zero + stat.n_part, dc_dt.begin() + HSO3 * stat.n_part, 
        rhs<HSO3, decltype(*this)>(*this)
      );
      thrust::transform(zero, zero + stat.n_part, dc_dt.begin() + SO3 * stat.n_part, 
        rhs<SO3, decltype(*this)>(*this)
      );
    }
  };
}
