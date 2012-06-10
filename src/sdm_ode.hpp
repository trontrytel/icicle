/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date March-June 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#pragma once

#include "sdm_base.hpp"

#if defined(USE_THRUST) && defined(USE_BOOST_ODEINT)
#  include <boost/numeric/odeint.hpp>
#  include <boost/numeric/odeint/external/thrust/thrust_algebra.hpp>
#  include <boost/numeric/odeint/external/thrust/thrust_operations.hpp>
#  include <boost/numeric/odeint/external/thrust/thrust_resize.hpp>
namespace odeint = boost::numeric::odeint;

namespace sdm
{

/*
  On May 9, 2012, at 7:44 PM, Karsten Ahnert wrote:
  > ... unfortunately the Rosenbrock method cannot be used with any other state type than ublas.matrix.
  > ... I think, the best steppers for stiff systems and thrust are the
  > runge_kutta_fehlberg78 or the bulirsch_stoer with a very high order. But
  > should benchmark both steppers and choose the faster one.
*/

  typedef odeint::euler<
    thrust::device_vector<thrust_real_t>, // state type
    thrust_real_t, // value_type
    thrust::device_vector<thrust_real_t>, // deriv type
    thrust_real_t, // time type
    odeint::thrust_algebra, 
    odeint::thrust_operations
  > algo_euler;

  typedef odeint::runge_kutta4<
    thrust::device_vector<thrust_real_t>, // state type
    thrust_real_t, // value_type
    thrust::device_vector<thrust_real_t>, // deriv type
    thrust_real_t, // time type
    odeint::thrust_algebra, 
    odeint::thrust_operations
  > algo_rk4;

  // nested class: 
  template <typename real_t>
  class ode 
  {
    public: virtual void init() {}
    public: virtual void advance(
      thrust::device_vector<thrust_real_t> &x, 
      const quantity<si::time, real_t> &dt
    ) = 0;

    // TODO: a better name needed...
    public: virtual thrust_real_t transform(const thrust_real_t &x)
    {
      return x; // i.e. identity 
    }
  };

  // nested class: 
  template <typename real_t, class algo> 
  class ode_algo : public ode<real_t>
  {
    private: algo stepper;

    // pure virtual method
    public: virtual void operator()(
      const thrust::device_vector<thrust_real_t> &xy, 
      thrust::device_vector<thrust_real_t> &dxy_dt, 
      const thrust_real_t
    ) = 0;

    private: bool inited = false;
    protected: virtual void init() {}
 
    public: void advance(
      thrust::device_vector<thrust_real_t> &x, 
      const quantity<si::time, real_t> &dt
    )
    {
      // intended for calculations that cannot be placed in the constructor 
      // since at that time the initial values are not loaded yet to psi
      if (!inited) 
      {
        init();
        inited = true;
      }
      stepper.do_step(boost::ref(*this), x, 0, dt / si::seconds);
    }
  };

}
#endif
