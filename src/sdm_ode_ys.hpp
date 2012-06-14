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
  template <typename real_t, class algo>
  class ode_ys : public ode_algo<real_t, algo>
  { 
    // nested functor
    class term_vel
    {
      private: const stat_t<real_t> &stat;
      private: const sdm::ode<real_t> &F_xi;
      public: term_vel(
        const stat_t<real_t> &stat, 
        const sdm::ode<real_t> &F_xi
      ) : stat(stat), F_xi(F_xi) {}
      public: real_t operator()(const thrust_size_t id)
      {
        return - phc::vt<real_t>(
          real_t(F_xi.transform(stat.xi[id])) * si::metres,
          real_t(300) * si::kelvins,                   // TODO !!!
          real_t(1) * si::kilograms / si::cubic_metres // TODO !!!
        ) * si::seconds / si::metres;
      }
    };

    // private fields
    private: const stat_t<real_t> &stat;
    private: const envi_t<real_t> &envi;
    private: const sdm::ode<real_t> &F_xi;
    private: const thrust::counting_iterator<thrust_size_t> iter;
 
    // ctor 
    public: ode_ys(
      const stat_t<real_t> &stat, 
      const envi_t<real_t> &envi,
      const sdm::ode<real_t> &F_xi
    ) 
      : stat(stat), envi(envi), F_xi(F_xi), iter(0)
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
        term_vel(stat, F_xi)
      );
    }
  };
}
