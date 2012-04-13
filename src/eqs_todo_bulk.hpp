/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date April 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef EQS_TODO_BULK_HPP
#  define EQS_TODO_BULK_HPP

#  include "eqs_todo.hpp"
#  if defined(USE_BOOST_ODEINT)
#    include <boost/numeric/odeint.hpp>
namespace odeint = boost::numeric::odeint;
#  endif

template <typename real_t>
class eqs_todo_bulk : public eqs_todo<real_t> 
{
  // ctor
  public: eqs_todo_bulk(const grd<real_t> &grid)
    : eqs_todo<real_t>(grid, true)
  {}

  // RHS of the ODE to be solved
  private: class rhs
  { 
    private: real_t const *rhod;
    public: void init(const real_t *ptr, const real_t rhod_th, const real_t rhod_rv) 
    { 
      rhod = ptr; 
      update(rhod_th, rhod_rv);
    }

    public: quantity<phc::mixing_ratio, real_t> r;
    public: quantity<si::pressure, real_t> p;
    public: quantity<si::temperature, real_t> T;
    private: void update(const real_t rhod_th, const real_t rhod_rv)
    {
      r = rhod_rv / *rhod; 

      p = phc::p_1000<real_t>() * real_t(pow(
        (rhod_th * phc::R_d<real_t>() * si::kelvins * si::kilograms / si::cubic_metres) 
          / phc::p_1000<real_t>() * (real_t(1) + r / phc::eps<real_t>()),
        1 / (1 - phc::R_over_c_p(r))
      ));

      T = rhod_th / *rhod * si::kelvins * phc::exner<real_t>(p, r);
    }

    // F = d (rho_d * th) / d (rho_d * r) = rho_d ()
    public: void operator()(const real_t rhod_th, real_t &F, const real_t rhod_rv)
    {
      update(rhod_th, rhod_rv);
      F = - rhod_th * (
        phc::l_v<real_t>(T) / real_t(pow(real_t(1) + r, 2)) / phc::c_p(r) / T
        /* + ... TODO */
      );
    }
  };

  // the saturation adjustment (aka ,,bulk'' microphysics)
  public: void adjustments(
    int n, // TODO: mo¿e jednak bez n...
    const ptr_vector<mtx::arr<real_t>> &aux, 
    vector<ptr_vector<mtx::arr<real_t>>> &psi
  ) 
  {
#  if !defined(USE_BOOST_ODEINT)
    error_macro("eqs_todo_bulk requires icicle to be compiled with Boost.odeint");
#  else

    // state_type = value_type = deriv_type = time_type = real_t
    typedef odeint::runge_kutta4<real_t, real_t, real_t, real_t, 
      odeint::vector_space_algebra, odeint::default_operations, odeint::never_resizer> stepper;

    const mtx::arr<real_t>
      &rhod = aux[this->par.idx_rhod];
    mtx::arr<real_t>
      &rhod_rv = psi[this->par.idx_rhod_rv][n],
      &rhod_rl = psi[this->par.idx_rhod_rl][n],
      &rhod_th = psi[this->par.idx_rhod_th][n];

    stepper S; // TODO: in ctor (but thread safety!)
    //observer O(rhod_th, rhod_rv, rhod_rl);
    rhs F;

    for (int k = rhod.lbound(mtx::k); k <= rhod.ubound(mtx::k); ++k)
      for (int j = rhod.lbound(mtx::j); j <= rhod.ubound(mtx::j); ++j)
        for (int i = rhod.lbound(mtx::i); i <= rhod.ubound(mtx::i); ++i)
        {

          F.init(&rhod(i,j,k), rhod_th(i,j,k), rhod_rv(i,j,k));
          real_t eps = -.00001; // TODO!!!

          while (
            //rhod_rl(i,j,k) > eps ||  // TODO: evaporation
            rhod_rv(i,j,k) - eps > rhod(i,j,k) * phc::r_vs<real_t>(F.T, F.p)
          ) 
          {
            S.do_step(boost::ref(F), rhod_th(i,j,k), rhod_rv(i,j,k), eps);
            rhod_rl(i,j,k) -= eps;
            rhod_rv(i,j,k) += eps;
          }
        }
#  endif
  }
};
#endif
