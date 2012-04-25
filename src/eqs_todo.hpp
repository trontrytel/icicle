/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date March 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef EQS_TODO_HPP
#  define EQS_TODO_HPP

#  include "cmn.hpp"
#  include "eqs.hpp"
#  include "rhs_explicit.hpp"
#  include "grd.hpp"
#  include "phc.hpp"

template <typename real_t>
class eqs_todo : public eqs<real_t> 
{
  // nested class (well... struct)
  protected: struct params
  {
    quantity<si::mass_density, real_t> rho_unit;
    int idx_rhod, idx_rhod_rv, idx_rhod_th, idx_rhod_rl, idx_rhod_rr; // auxiliary variable indices
  };

  // nested class TODO this is only needed for bulk
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
 
      Ra(Ra.ijk) += sign*(max( 0., .001 * 
                    (rhod_rl(rhod.ijk) / rhod(rhod.ijk) - .0005))); //should be .001
      cout << "Ra min: " << min(Ra(Ra.ijk)) <<" Ra max: " << max(Ra(Ra.ijk)) << endl; 
      }
  };

  // nested class TODO this is only needed for bulk
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
      const quantity<si::time, real_t>
      )
      { 
      const mtx::arr<real_t>
        &rhod = aux[par.idx_rhod],
        &rhod_rl = (*psi[par.idx_rhod_rl]),
        &rhod_rr = (*psi[par.idx_rhod_rr]);
      
      assert(min(rhod_rr) >=0);
      assert(min(rhod_rl) >=0);

      Rc(Rc.ijk) += sign * 2.2 * rhod_rl(rhod.ijk) * 
                    ( pow(rhod_rr(rhod.ijk)/rhod(rhod.ijk), .875) );
      cout << "Rc min: " << min(Rc(Rc.ijk)) <<" Rc max: "<< max(Rc(Rc.ijk)) << endl;
      }
  };

  // nested class TODO this is only needed for bulk
/*  private: class rain_evaporation : public rhs_explicit<real_t>
  {
    private: struct params &par;
    private: real_t sign ;
    
    //ctor
    public: rain_evaporation(real_t sign, struct params &par)
    : sign(sign), par(par)
    {}

    public: void explicit_part(
      mtx::arr<real_t> &Rc,
      const ptr_vector<mtx::arr<real_t>> &aux,
      const mtx::arr<real_t> * const * const psi,
      const quantity<si::time, real_t>
      )
      { 
      const mtx::arr<real_t>
        &rhod = aux[par.idx_rhod],
        &rhod_rl = (*psi[par.idx_rhod_rl]),
        &rhod_rr = (*psi[par.idx_rhod_rr]);

      Re(Re.ijk) += 0;
      }
  };
*/
 
  // nested class TODO this is only needed for bulk
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
     cout << "min rhod_rr  "<<min(rhod_rr)<<endl;     
     cout << "max rhod_rr  "<<max(rhod_rr)<<endl;     
 
      //TODO do something about fabs() for rhod_rr
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
/*
              ((rhod_rr(i,j+1,k) * 36.34 * pow(rhod_rr(i,j+1,k)*1e-3,.1346) 
              * pow(rhod(i,rhod.lbound(mtx::j),k)/rhod(i,j+1,k),.5) )
              -
              (rhod_rr(i,j,k) * 36.34 * pow(rhod_rr(i,j,k)*1e-3,.1346) 
              * pow(rhod(i,rhod.lbound(mtx::j),k)/rhod(i,j,k),.5) ))
              /real_t(grid.dy() / si::metres ) ;
*/
          }
          Rt(i, rhod.lbound(mtx::j), k) -= real_t(
            rhod_rr(i, rhod.lbound(mtx::j)+1,k)-rhod_rr(i, rhod.lbound(mtx::j),k)
            ) * v_iph / si::metres_per_second / real_t(grid.dy() / si::metres)  ;
        }
      }
      cout << "Rt min: " << min(Rt(Rt.ijk)) << " Rt max: " << max(Rt(Rt.ijk)) << endl;
    }
  };

  // private field
  protected: params par;

  // ctor
  public: eqs_todo(const grd<real_t> &grid, bool bulk)
  {
    par.rho_unit = 1 * si::kilograms / si::cubic_metres;

    // auxiliary variable
    this->aux.push_back(new struct eqs<real_t>::axv({
      "rhod", "dry air density", this->quan2str(par.rho_unit),
      vector<int>({0, 0, 0})
    }));
    par.idx_rhod = this->aux.size() - 1;

    // eg. eqn 1b in Szumowski, Grabowski & Ochs 1998, Atmos. Res.
    // cf. eqn 3.55 in the Jacobson's Fundamentals of Atmos. Modelling (2nd edn, 2005)
    this->sys.push_back(new struct eqs<real_t>::gte({
      "rhod_rv", "dry air density times water vapour mixing ratio (i.e. water vapour mass density or absolute humidity)",
      this->quan2str(par.rho_unit), 
      typename eqs<real_t>::positive_definite(true),
    }));
    par.idx_rhod_rv = this->sys.size() - 1;

    // only for bulk model (i.e. not for super droplets)
    //TODO asserts for rl and rr ge 0
    if (bulk) // TODO: make it an aux var for SDM?
    {
      this->sys.push_back(new struct eqs<real_t>::gte({
        "rhod_rl", "dry air density times liquid water mixing ratio (i.e. liquid water mass density)",
        this->quan2str(par.rho_unit), 
        typename eqs<real_t>::positive_definite(true),
      }));
      this->sys.back().rhs_terms.push_back(new collection(-1, par));
      this->sys.back().rhs_terms.push_back(new autoconversion(-1, par));
      par.idx_rhod_rl = this->sys.size() - 1;

      this->sys.push_back(new struct eqs<real_t>::gte({
        "rhod_rr", "dry air density times rain water mixing ratio (i.e. rain water mass density)",
        this->quan2str(par.rho_unit), 
        typename eqs<real_t>::positive_definite(true),
      }));
      //TODO should not be calculated twice
      this->sys.back().rhs_terms.push_back(new terminal_velocity(par, grid));
      this->sys.back().rhs_terms.push_back(new collection(+1, par));
      this->sys.back().rhs_terms.push_back(new autoconversion(+1, par));
      par.idx_rhod_rr = this->sys.size() - 1;
    }

    // eg. eqn 1a in Szumowski, Grabowski & Ochs 1998, Atmos. Res.
    // cf. eqn 3.68 in the Jacobson's Fundamentals of Atmos. Modelling (2nd edn, 2005)
    this->sys.push_back(new struct eqs<real_t>::gte({
      "rhod_th", "dry air density times potential temperature (i.e. energy mass density over specific heat capacity)",
      this->quan2str(par.rho_unit * si::kelvins), 
      typename eqs<real_t>::positive_definite(true),
    }));
    par.idx_rhod_th = this->sys.size() - 1;
  }
};
#endif
