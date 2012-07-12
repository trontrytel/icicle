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
    // nested functor 
    class rhs : public xi
    {
      private: const stat_t<real_t> &stat;
      private: const envi_t<real_t> &envi;
      public: rhs(
        const stat_t<real_t> &stat,
        const envi_t<real_t> &envi
      ) : stat(stat), envi(envi) {}

      // overloaded operator invoked by thrust
      public: real_t operator()(const thrust_size_t id)
      {
cerr << "CHEMISTRY!!!" << endl;
        //// Aniu! to moze sie przydac:
        // thrust_size_t ij = stat.ij[id];
        // this->rw_of_xi(stat.xi[id]) 
        // envi.T[ij] 
        // envi.rhod[ij]
        return 0; // ta funkcja ma zwrocic wartosci "prawych stron" ODE
      }
    };

    // private fields
    private: const stat_t<real_t> &stat;
    private: const envi_t<real_t> &envi;
    private: const thrust::counting_iterator<thrust_size_t> iter;
 
    // ctor 
    public: ode_chem(
      const stat_t<real_t> &stat, 
      const envi_t<real_t> &envi
    ) 
      : stat(stat), envi(envi), iter(0)
    {}

    // overloaded () operator invoked by odeint
    public: void operator()(
      const thrust::device_vector<real_t>&, 
      thrust::device_vector<real_t> &dc_dt, 
      const real_t
    )
    {
      thrust::transform(
        iter, iter + 4 * stat.n_part, // TODO-AJ: 4 - od czterech stezen!
        dc_dt.begin(),  
        rhs(stat, envi)
      );
    }
  };
}
