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

    public: quantity<phc::mixing_ratio, real_t> r, rs; 
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

      rs = phc::r_vs<real_t>(T, p);
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
    vector<ptr_vector<mtx::arr<real_t>>> &psi,
    const quantity<si::time, real_t> dt
  ) 
  {
    bool  // TODO: take values from options
      opt_rainevap = false, 
      opt_cloudevap = true;
//      opt_autoconversion = true,
//      opt_sedimentation = true;

    const mtx::arr<real_t>
      &rhod = aux[this->par.idx_rhod];
    mtx::arr<real_t>
      &rhod_th = psi[this->par.idx_rhod_th][n],
      &rhod_rv = psi[this->par.idx_rhod_rv][n],
      &rhod_rl = psi[this->par.idx_rhod_rl][n],
      &rhod_rr = psi[this->par.idx_rhod_rr][n];

//    if (opt_collection) collection + limiter
//    if (opt_autoconversion) autoconversion(rhod, rhod_rl, rhod_rr, dt);
//    if (opt_sedimentation) sedimentation(rhod, rhod_rr, dt);
    condevap(n, aux, psi, dt, opt_rainevap, opt_cloudevap);
  }
 
/*   private: void sedimentation(
    const mtx::arr<real_t> &rhod,
    mtx::arr<real_t> &rhod_rr,
    const quantity<si::time, real_t> dt
  )
  {
    // TODO! (needs temp var!)
  }

  private: void autoconversion(
    const mtx::arr<real_t> &rhod,
    mtx::arr<real_t> &rhod_rl,
    mtx::arr<real_t> &rhod_rr,
    const quantity<si::time, real_t> dt
  )
  {
    assert(min(rhod_rl) >=0);
    assert(min(rhod_rr) >=0);
 
    rhod_rl(rhod.ijk) -= (dt / si::seconds) * (max( 0., .001 * (rhod_rl(rhod.ijk) / rhod(rhod.ijk) - .0005))); //should be .001
    rhod_rr(rhod.ijk) += (dt / si::seconds) * (max( 0., .001 * (rhod_rl(rhod.ijk) / rhod(rhod.ijk) - .0005))); //should be .001
    //cout << "Ra min: " << min(Ra(Ra.ijk)) <<" Ra max: " << max(Ra(Ra.ijk)) << endl; 
  }
*/ 
  private: void condevap(
    int n, // TODO: mo¿e jednak bez n...
    const ptr_vector<mtx::arr<real_t>> &aux, 
    vector<ptr_vector<mtx::arr<real_t>>> &psi,
    const quantity<si::time, real_t> dt,
    bool opt_rainevap, 
    bool opt_cloudevap
  )   
  {
#  if !defined(USE_BOOST_ODEINT)
    error_macro("eqs_todo_bulk requires icicle to be compiled with Boost.odeint");
#  else

    const mtx::arr<real_t>
      &rhod = aux[this->par.idx_rhod];
    mtx::arr<real_t>
      &rhod_rl = psi[this->par.idx_rhod_rl][n],
      &rhod_rr = psi[this->par.idx_rhod_rr][n],
      &rhod_rv = psi[this->par.idx_rhod_rv][n],
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
          real_t drho_rr_max = 0; // TODO: quantity<si::mass_density
          if (F.rs > F.r && rhod_rr(i,j,k) > 0 && opt_rainevap) 
            drho_rr_max = (dt / si::seconds) * (1 - F.r / F.rs) * (1.6 + 124.9 * pow(1e-3 * rhod_rr(i,j,k), .2046)) * pow(1e-3 * rhod_rr(i,j,k), .525) / 
              (5.4e2 + 2.55e5 * (1. / (F.p / si::pascals) / F.rs));
          bool incloud;

          // TODO: rethink and document 2*rho_eps!!!
          while ( 
            // condensation of cloud water if supersaturated
            (vapour_excess = rhod_rv(i,j,k) - rhod(i,j,k) * F.rs) > rho_eps 
   //         || // or ...
   //         (vapour_excess < -rho_eps && ( // if subsaturated
   //           (incloud = (rhod_rl(i,j,k) > 0)) // cloud evaporation if in cloud
   //           || (opt_rainevap && rhod_rr(i,j,k) > 0) // or rain evaportation if in rain shaft (and out-of-cloud)
   //         )) 
          ) 
          {
            real_t drho_rv = - copysign(.5 * rho_eps, vapour_excess);
            drho_rv = (vapour_excess > 0 || incloud)
              ?  std::min(rhod_rl(i,j,k), drho_rv)
              :  std::min(drho_rr_max, std::min(rhod_rr(i,j,k), drho_rv)); // preventing negative mixing ratios
            assert(drho_rv != 0); // otherwise it should not pass the while condition!

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

            // latent heat source/sink due to evaporation/condensation
            rhod_th(i,j,k) = tmp / (si::kilograms / si::cubic_metres * si::kelvins); 

            // updating rhod_rv
            rhod_rv(i,j,k) += drho_rv;
            assert(rhod_rv(i,j,k) >= 0);
            assert(isfinite(rhod_rv(i,j,k)));
            
            if (vapour_excess > 0 || incloud) 
            {
              rhod_rl(i,j,k) -= drho_rv; // cloud water 
              assert(rhod_rl(i,j,k) >= 0);
              assert(isfinite(rhod_rl(i,j,k)));
            }
            else // or rain water
            {
              assert(opt_rainevap);
              rhod_rr(i,j,k) -= drho_rv;
              assert(rhod_rr(i,j,k) >= 0);
              assert(isfinite(rhod_rr(i,j,k)));
              if ((drho_rr_max -= drho_rv) == 0) break; // but not more than Kessler allows
            }
          }
          // hopefully true for RK4
          assert(F.r == real_t(rhod_rv(i,j,k) / rhod(i,j,k)));
          // double-checking....
          assert(rhod_rl(i,j,k) >= 0);
          assert(rhod_rv(i,j,k) >= 0);
          assert(rhod_rr(i,j,k) >= 0);
        }
#  endif
  }
};
#endif
