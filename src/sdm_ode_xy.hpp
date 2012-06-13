/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date March-June 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#pragma once
#include "sdm_ode.hpp"

namespace sdm
{
  template <typename real_t, class algo>
  class ode_xy : public ode_algo<real_t, algo>
  { 
    // nested functor 
    private: 
    template <int di, int dj>
    class interpol
    {
      // private fields
      private: const int n;
      private: const thrust::device_vector<thrust_real_t> &vel;
      private: const stat_t<real_t> &stat;

      // ctor
      public: interpol(
        const int n,
        const thrust::device_vector<thrust_real_t> &vel,
        const stat_t<real_t> &stat
      ) : n(n), vel(vel), stat(stat) {}

      //  overloaded () operator invoked by thrust::transform()
      public: thrust_real_t operator()(thrust_size_t id)
      {
        // TODO weighting by position!
        return .5 * (
          vel[(stat.i[id]     ) + (stat.j[id]     ) * n] + 
          vel[(stat.i[id] + di) + (stat.j[id] + dj) * n]
        );
      }
    };

    // private fields
    private: const stat_t<real_t> &stat;
    private: const envi_t<real_t> &envi;
 
    // ctor 
    public: ode_xy(const stat_t<real_t> &stat, const envi_t<real_t> &envi) : stat(stat), envi(envi) {}

    // overloaded () operator invoked by odeint
    public: void operator()(
      const thrust::device_vector<thrust_real_t>&, 
      thrust::device_vector<thrust_real_t> &dxy_dt, 
      const thrust_real_t
    )
    {
      // TODO use positions to interpolate velocities! (as an option?)
      thrust::counting_iterator<thrust_size_t> iter(0);
      thrust::transform(
        iter, iter + stat.n_part, 
        dxy_dt.begin(), 
        interpol<1,0>(envi.vx_nx, envi.vx, stat)
      );
      thrust::transform(
        iter, iter + stat.n_part, 
        dxy_dt.begin() + stat.n_part, 
        interpol<0,1>(envi.vy_nx, envi.vy, stat)
      );
    }
  };
}
