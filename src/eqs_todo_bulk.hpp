/** @file
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date April 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef EQS_TODO_BULK_HPP
#  define EQS_TODO_BULK_HPP

#  include "eqs_todo.hpp"
#  include "phc_kessler.hpp"
#  include "phc_theta.hpp"

template <typename real_t>
class eqs_todo_bulk : public eqs_todo<real_t> 
{
  // nested class (well... struct)
  protected: struct params : eqs_todo<real_t>::params
  {
    int idx_rhod_rl, idx_rhod_rr; // advected variables indices
  };
  protected: params par;

  // a container for storing options (i.e. which processes ar on/off)
  public: enum processes {cond, cevp, conv, clct, sedi, revp};
  private: map<enum processes, bool> opts;

  // nested class
  private: class autoconversion : public rhs_explicit<real_t>
  {
    private: struct params &par;
    private: real_t sign ;

    //ctor
    public: autoconversion(real_t sign, struct params &par)
    : sign(sign), par(par)
    {}

    public: void explicit_part(
      mtx::arr<real_t> &Ra,
      const ptr_vector<mtx::arr<real_t>> &aux,
      const mtx::arr<real_t> * const * const psi,
      const quantity<si::time, real_t>
      )
    {
      const mtx::arr<real_t>
        &rhod = aux[par.idx_rhod],
        &rhod_rl = (*psi[par.idx_rhod_rl]);

      assert(min(rhod_rl) >=0);

      // TODO: threshold value as an option!
      // TODO: move the formulat to phc
      Ra(Ra.ijk) += sign*(max( 0., .001 * rhod(rhod.ijk) * (rhod_rl(rhod.ijk) / rhod(rhod.ijk) - .0005))); //should be .001
    }
  };

  // nested class
  private: class collection : public rhs_explicit<real_t>
  {
    private: struct params &par;
    private: real_t sign ;

    //ctor
    public: collection(real_t sign, struct params &par)
    : sign(sign), par(par)
    {}

    public: void explicit_part(
      mtx::arr<real_t> &Rc,
      const ptr_vector<mtx::arr<real_t>> &aux,
      const mtx::arr<real_t> * const * const psi,
      const quantity<si::time, real_t> dt
    )
    {
      const mtx::arr<real_t>
        &rhod = aux[par.idx_rhod],
        &rhod_rl = (*psi[par.idx_rhod_rl]),
        &rhod_rr = (*psi[par.idx_rhod_rr]);

      assert(min(rhod_rr) >=0);
      assert(min(rhod_rl) >=0);

      // TODO: move the formulat to phc
      Rc(Rc.ijk) += sign * 2.2 * rhod_rl(rhod.ijk) * ( pow(rhod_rr(rhod.ijk)/rhod(rhod.ijk), .875) );
    }
  };

  // nested class
  private: class terminal_velocity : public rhs_explicit<real_t>
  {
    private: struct params &par;
    private: const grd<real_t> &grid;

    //ctor
    public: terminal_velocity(struct params &par, const grd<real_t> &grid)
    : par(par), grid(grid)
    {}

    public: void explicit_part(
      mtx::arr<real_t> &Rt,
      const ptr_vector<mtx::arr<real_t>> &aux,
      const mtx::arr<real_t> * const * const psi,
      const quantity<si::time, real_t>
    )
    {
      const mtx::arr<real_t>
        &rhod = aux[par.idx_rhod],
        &rhod_rr = (*psi[par.idx_rhod_rr]);

      assert(min(rhod_rr) >=0);

      for (int k = rhod.lbound(mtx::k); k <= rhod.ubound(mtx::k); ++k)
      {
        for (int i = rhod.lbound(mtx::i); i <= rhod.ubound(mtx::i); ++i)
        {
          quantity<si::velocity, real_t> v_iph = 0 * si::metres_per_second;
          for (int j = rhod.ubound(mtx::j); j > rhod.lbound(mtx::j); --j)
          {  //upstream?
            quantity<si::velocity, real_t> v_imh = real_t(-.5) * ( // -.5 : averaging + axis orientation
              phc::vterm<real_t>(
                rhod_rr(i,j-1,k) * si::kilograms / si::cubic_metres,
                rhod(i,j-1,k) * si::kilograms / si::cubic_metres,
                rhod(i,rhod.lbound(mtx::j),k) * si::kilograms / si::cubic_metres
              ) +
              phc::vterm<real_t>(
                rhod_rr(i,j,k) * si::kilograms / si::cubic_metres,
                rhod(i,j,k) * si::kilograms / si::cubic_metres,
                rhod(i,rhod.lbound(mtx::j),k) * si::kilograms / si::cubic_metres
              ));
#define upstream_F(a,b,c) ((b)*(c)) // assumes c is negative!
            Rt(i,j,k) -= real_t(
              upstream_F(rhod_rr(i,j,  k), rhod_rr(i,j+1,k), v_iph / si::metres_per_second) -
              upstream_F(rhod_rr(i,j-1,k), rhod_rr(i,j,  k), v_imh / si::metres_per_second)
            ) / real_t(grid.dy() / si::metres);
            v_iph = v_imh;
          }
          Rt(i, rhod.lbound(mtx::j), k) -= real_t(
            rhod_rr(i, rhod.lbound(mtx::j)+1,k)-rhod_rr(i, rhod.lbound(mtx::j),k)
            ) * v_iph / si::metres_per_second / real_t(grid.dy() / si::metres)  ;
// TODO: record accumulated rainfall
        }
      }
    }
  };

  // ctor
  public: eqs_todo_bulk(const grd<real_t> &grid, map<enum processes, bool> opts)
    : eqs_todo<real_t>(grid, &par), opts(opts)
  {
    if (opts[revp] && !opts[cevp])
      error_macro("rain evaporation requires cloud evaporation to be turned on")

    this->sys.push_back(new struct eqs<real_t>::gte({
      "rhod_rl", "dry air density times liquid water mixing ratio (i.e. liquid water mass density)",
      this->quan2str(this->par.rho_unit),
      typename eqs<real_t>::positive_definite(true),
    }));
    if (opts[clct]) this->sys.back().rhs_terms.push_back(new collection(-1, par)); 
    if (opts[conv]) this->sys.back().rhs_terms.push_back(new autoconversion(-1, par)); 
    par.idx_rhod_rl = this->sys.size() - 1;

    this->sys.push_back(new struct eqs<real_t>::gte({
      "rhod_rr", "dry air density times rain water mixing ratio (i.e. rain water mass density)",
      this->quan2str(par.rho_unit),
      typename eqs<real_t>::positive_definite(true),
    }));
    //TODO should not be calculated twice
    if (opts[clct]) this->sys.back().rhs_terms.push_back(new collection(+1, par));
    if (opts[conv]) this->sys.back().rhs_terms.push_back(new autoconversion(+1, par));
    if (opts[sedi]) this->sys.back().rhs_terms.push_back(new terminal_velocity(par, grid)); 
    par.idx_rhod_rr = this->sys.size() - 1;
  }

  // the saturation adjustment (aka ,,bulk'' microphysics)
  public: void adjustments(
    int n, // TODO: mo¿e jednak bez n...
    vector<ptr_vector<mtx::arr<real_t>>> &psi,
    ptr_vector<mtx::arr<real_t>> &aux, 
    const ptr_vector<mtx::arr<real_t>> C,
    const quantity<si::time, real_t> dt
  ) 
  {
    const mtx::arr<real_t>
      &rhod = aux[this->par.idx_rhod];
    mtx::arr<real_t>
      &rhod_th = psi[this->par.idx_rhod_th][n],
      &rhod_rv = psi[this->par.idx_rhod_rv][n],
      &rhod_rl = psi[this->par.idx_rhod_rl][n],
      &rhod_rr = psi[this->par.idx_rhod_rr][n];

    condevap(rhod, rhod_th, rhod_rv, rhod_rl, rhod_rr, dt);
  }
 
  private: virtual void condevap(
    const mtx::arr<real_t> &rhod,
    mtx::arr<real_t> &rhod_th,
    mtx::arr<real_t> &rhod_rv,
    mtx::arr<real_t> &rhod_rl,
    mtx::arr<real_t> &rhod_rr,
    const quantity<si::time, real_t> dt
  ) = 0;   
};
#endif
