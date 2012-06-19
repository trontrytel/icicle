/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date June 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#pragma once
#include "sdm_ode.hpp"
#include "phc_terminal_vel.hpp"

namespace sdm
{
  template <typename real_t, class algo, class xi>
  class ode_ys : public ode_algo<real_t, algo>
  { 
    // nested functor
    class term_vel : public xi
    {
      private: const stat_t<real_t> &stat;
      private: const envi_t<real_t> &envi;
      public: term_vel(
        const stat_t<real_t> &stat,
        const envi_t<real_t> &envi
      ) : stat(stat), envi(envi) {}
      public: real_t operator()(const thrust_size_t id)
      {
        return - phc::vt<real_t>(
          this->rw_of_xi(stat.xi[id]) * si::metres,
          envi.T[stat.ij[id]] * si::kelvins,
          envi.rhod[stat.ij[id]] * si::kilograms / si::cubic_metres // that's the dry air density
        ) * si::seconds / si::metres;
      }
    };

    // private fields
    private: const stat_t<real_t> &stat;
    private: const envi_t<real_t> &envi;
    private: const thrust::counting_iterator<thrust_size_t> iter;
 
    // ctor 
    public: ode_ys(
      const stat_t<real_t> &stat, 
      const envi_t<real_t> &envi
    ) 
      : stat(stat), envi(envi), iter(0)
    {}

    // overloaded () operator invoked by odeint
    public: void operator()(
      const thrust::device_vector<real_t>&, 
      thrust::device_vector<real_t> &dxy_dt, 
      const real_t
    )
    {
      thrust::transform(
        iter, iter + stat.n_part, 
        dxy_dt.begin() + stat.n_part, 
        term_vel(stat, envi)
      );
    }
  };
}
