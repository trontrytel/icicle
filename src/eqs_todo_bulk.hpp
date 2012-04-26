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
    private: quantity<si::mass_density, real_t> rhod;
    public: void init(
      const quantity<si::mass_density, real_t> _rhod, 
      const quantity<multiply_typeof_helper<si::mass_density, si::temperature>::type, real_t> rhod_th, 
      const quantity<si::mass_density, real_t> rhod_rv
    ) 
    { 
      rhod = _rhod; 
      update(rhod_th, rhod_rv);
    }

    public: quantity<phc::mixing_ratio, real_t> r; // TODO: nie! r jest przecie� zmienn�!!!
    public: quantity<si::pressure, real_t> p;
    public: quantity<si::temperature, real_t> T;
    private: void update(
      const quantity<multiply_typeof_helper<si::mass_density, si::temperature>::type, real_t> rhod_th, 
      const quantity<si::mass_density, real_t> rhod_rv
    )
    {
      r = rhod_rv / rhod; 

      p = phc::p_1000<real_t>() * real_t(pow(
        (rhod_th * phc::R_d<real_t>()) 
          / phc::p_1000<real_t>() * (real_t(1) + r / phc::eps<real_t>()),
        1 / (1 - phc::R_over_c_p(r))
      ));

      T = rhod_th / rhod * phc::exner<real_t>(p, r);
    }

    // F = d (rho_d * th) / d (rho_d * r) = rho_d ()
    public: void operator()(
      const quantity<multiply_typeof_helper<si::mass_density, si::temperature>::type, real_t> rhod_th, 
      quantity<si::temperature, real_t> &F, 
      const quantity<si::mass_density, real_t> rhod_rv
    )
    {
      update(rhod_th, rhod_rv);
      F = - rhod_th / rhod * (
        phc::l_v<real_t>(T) / real_t(pow(1 + r, 2)) / phc::c_p(r) / T // the 'liquid water' term
        //+ 
        //log(p/phc::p_1000<real_t>()) * phc::R_d_over_c_pd<real_t>() * (1/phc::eps<real_t>() - 1/phc::ups<real_t>()) * pow(1+r/phc::ups<real_t>(),-2) // the 'virtual' term
      );
//cerr << (
//        log(p/phc::p_1000<real_t>()) * phc::R_d_over_c_pd<real_t>() * (1/phc::eps<real_t>() - 1/phc::ups<real_t>()) * pow(1+r/phc::ups<real_t>(),-2) // the 'virtual' term
//) / (
//        phc::l_v<real_t>(T) / pow(1 + r, 2) / phc::c_p(r) / T // the 'liquid water' term
//) << endl;
    }
  };

  //typedef odeint::euler< // TODO: opcja?
  typedef odeint::runge_kutta4<
    quantity<multiply_typeof_helper<si::mass_density, si::temperature>::type, real_t>, // state_type
    real_t, // value_type
    quantity<si::temperature, real_t>, // deriv_type
    quantity<si::mass_density, real_t>, // time_type
    odeint::vector_space_algebra, 
    odeint::default_operations, 
    odeint::never_resizer
  > stepper;

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

    const mtx::arr<real_t>
      &rhod = aux[this->par.idx_rhod];
    mtx::arr<real_t>
      &rhod_rv = psi[this->par.idx_rhod_rv][n],
      &rhod_rl = psi[this->par.idx_rhod_rl][n],
      &rhod_rr = psi[this->par.idx_rhod_rr][n],
      &rhod_th = psi[this->par.idx_rhod_th][n];

    stepper S; // TODO: in ctor (but thread safety!)
    rhs F;

    for (int k = rhod.lbound(mtx::k); k <= rhod.ubound(mtx::k); ++k)
      for (int j = rhod.lbound(mtx::j); j <= rhod.ubound(mtx::j); ++j)
        for (int i = rhod.lbound(mtx::i); i <= rhod.ubound(mtx::i); ++i)
        {

          F.init(
            rhod(i,j,k)    * si::kilograms / si::cubic_metres, 
            rhod_th(i,j,k) * si::kilograms / si::cubic_metres * si::kelvins, 
            rhod_rv(i,j,k) * si::kilograms / si::cubic_metres
          );
          real_t // TODO: quantity<si::mass_density
            rho_eps = .00002, // TODO: as an option?
            vapour_excess;
          real_t // TODO: quantity<si::mass_density
            drho_rr_max = 0.01; //TODO!
          bool incloud;

          // TODO: something more meaningfull than 2*rho_eps!!!
          while ( 
            // condensation of cloud water if supersaturated
            (vapour_excess = rhod_rv(i,j,k) - rhod(i,j,k) * phc::r_vs<real_t>(F.T, F.p)) > rho_eps 
            || // or ...
            (vapour_excess < -rho_eps && ( // if subsaturated
              (incloud = (rhod_rl(i,j,k) > 10 * mtx::eps<real_t>())) // cloud evaporation if in cloud
              || rhod_rr(i,j,k) > 10 * mtx::eps<real_t>() // or rain evaportation if in rain shaft (and out-of-cloud)
            ))  // TODO: brzydkie eps!, dzialalo dla >0 dla samej chmury!
          ) 
          {
            real_t drho_rv = - copysign(.5 * rho_eps, vapour_excess);
            drho_rv = (vapour_excess > 0 || incloud)
              ?  std::min(rhod_rl(i,j,k), drho_rv)
              :  std::min(drho_rr_max, std::min(rhod_rr(i,j,k), drho_rv)); // preventing negative mixing ratios
            //if (abs(drho_rv) < mtx::eps<real_t>()) break;// TODO: moze jako assert()?

cerr << "rhod_rl=" << rhod_rl(i,j,k) << ", rhod_rr=" << rhod_rr(i,j,k) << ", drho=" << drho_rv << endl;

            // theta is modified by do_step, and hence we cannot pass an expression and we need a temp. var.
            quantity<multiply_typeof_helper<si::mass_density, si::temperature>::type, real_t> 
              tmp = rhod_th(i,j,k) * si::kilograms / si::cubic_metres * si::kelvins;
            // integrating the First Law for moist air
            S.do_step(
              boost::ref(F), 
              tmp,
              rhod_rv(i,j,k) * si::kilograms / si::cubic_metres, 
              drho_rv        * si::kilograms / si::cubic_metres
            );
            // latent heat source/sink due to ...
            rhod_th(i,j,k) = tmp / (si::kilograms / si::cubic_metres * si::kelvins); 
            // ... condensation/evaporation of ...
            rhod_rv(i,j,k) += drho_rv;
            if (vapour_excess > 0 || incloud) 
              rhod_rl(i,j,k) -= drho_rv; // cloud water 
            else // or rain water
            {
              rhod_rr(i,j,k) -= drho_rv;
              if ((drho_rr_max -= drho_rv) == 0) break; // but not more than Kessler allows
            }
          }
          // hopefully true for RK4
          assert(F.r == real_t(rhod_rv(i,j,k) / rhod(i,j,k)));
        }
#  endif
  }
};
#endif
