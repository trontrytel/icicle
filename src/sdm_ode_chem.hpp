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
//#include "phc_henry...hpp" // TODO-AJ
//#include "phc_react...hpp" // TODO-AJ

namespace sdm
{
  template <typename real_t, class algo, class xi>
  class ode_chem : public ode_algo<real_t, algo>
  { 
    enum {SO2, HSO3, SO3, H2O2};
        //// Aniu! to moze sie przydac:
        // thrust_size_t ij = stat.ij[id];
        // this->rw_of_xi(stat.xi[id]) 
        // envi.T[ij] 
        // envi.rhod[ij]

    // nested functor 
    template <int chem>
    class rhs : public xi
    {
      private: const stat_t<real_t> &stat;
      private: const envi_t<real_t> &envi;
      public: rhs(const stat_t<real_t> &stat, const envi_t<real_t> &envi) 
        : stat(stat), envi(envi) {}

      // overloaded operator invoked by thrust
      public: real_t operator()(const thrust_size_t id)
      {
        switch (chem)
        {
          case SO2: return 0; //(c_equil - stat.c[id + n_part * SO2]) / dt; 
          case HSO3: return 0; 
          case SO3: return 0; 
          case H2O2: return 0; 
        }
      }
    };

    // private fields
    private: const stat_t<real_t> &stat;
    private: const envi_t<real_t> &envi;
    private: const thrust::counting_iterator<thrust_size_t> zero;
    private: quantity<si::time, real_t> dt;
 
    void init(const quantity<si::time, real_t> &dt_arg)
    {
      dt = dt_arg;
    }

    // ctor 
    public: ode_chem(
      const stat_t<real_t> &stat, 
      const envi_t<real_t> &envi,
      const real_t c_SO2,
      const real_t c_HSO3,
      const real_t c_SO3,
      const real_t c_H2O2
    ) 
      : stat(stat), envi(envi), zero(0)
    {}

    // overloaded () operator invoked by odeint
    public: void operator()(
      const thrust::device_vector<real_t> &, 
      thrust::device_vector<real_t> &dc_dt, 
      const real_t
    )
    {
      thrust::transform(zero, zero + stat.n_part, dc_dt.begin() + SO2  * stat.n_part, rhs<SO2>(stat, envi));
      thrust::transform(zero, zero + stat.n_part, dc_dt.begin() + HSO3 * stat.n_part, rhs<HSO3>(stat, envi));
      thrust::transform(zero, zero + stat.n_part, dc_dt.begin() + SO3  * stat.n_part, rhs<SO3>(stat, envi));
      thrust::transform(zero, zero + stat.n_part, dc_dt.begin() + H2O2 * stat.n_part, rhs<H2O2>(stat, envi));
      //thrust::transform(zero, zero + stat.n_part, dc_dt.begin() + 4 * stat.n_part, rhs_O3(stat, envi)); // TODO
    }
  };
}
