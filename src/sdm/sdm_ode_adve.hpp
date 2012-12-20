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
  class ode_adve : public ode_algo<real_t, algo>
  { 
    // nested functor 
    private: 
    template <int if_x, int if_y> // TODO: int -> bool?
    class interpol
    {
      // private fields
      private: const thrust_size_t n;
      private: const thrust::device_vector<real_t> &vel;
      private: const real_t dxy;
      private: const stat_t<real_t> &stat;

      // ctor
      public: interpol(
        const thrust_size_t n,
        const thrust::device_vector<real_t> &vel,
        const real_t &dxy,
        const stat_t<real_t> &stat
      ) : n(n), vel(vel), stat(stat), dxy(dxy) {}

      //  overloaded () operator invoked by thrust::transform()
      public: real_t operator()(thrust_size_t id)
      {
        // weighting by position (assumes Arakawa-C grid)
        real_t w = 
        (
          (
            if_x * *(stat.x_begin + id) +
            if_y * *(stat.y_begin + id)
          ) - (
            if_x * stat.i[id] + 
            if_y * stat.j[id]
          ) * dxy
        ) / (dxy);
        return 
          (1 - w) * vel[(stat.i[id]       ) + (stat.j[id]       ) * n] + 
             w    * vel[(stat.i[id] + if_x) + (stat.j[id] + if_y) * n]
        ;
      }
    };

    // private fields
    private: const stat_t<real_t> &stat;
    private: const envi_t<real_t> &envi;
 
    // ctor 
    public: ode_adve(const stat_t<real_t> &stat, const envi_t<real_t> &envi) : stat(stat), envi(envi) {}

    // overloaded () operator invoked by odeint
    public: void operator()(
      const thrust::device_vector<real_t>&, 
      thrust::device_vector<real_t> &dxy_dt, 
      const real_t t
    )
    {
      // TODO: adjust() method to re-interpolate velocities after each sub-step (as vapour in cond)?
      thrust::counting_iterator<thrust_size_t> iter(0);
      thrust::transform(
        iter, iter + stat.n_part, 
        dxy_dt.begin(), 
        interpol<1,0>(envi.vx_nx, envi.vx, envi.dx, stat)
      );
      thrust::transform(
        iter, iter + stat.n_part, 
        dxy_dt.begin() + stat.n_part, 
        interpol<0,1>(envi.vy_nx, envi.vy, envi.dy, stat)
      );
    }
  };
}
