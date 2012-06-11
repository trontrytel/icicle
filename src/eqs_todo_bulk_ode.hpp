/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date April 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#pragma once
#include "eqs_todo_bulk.hpp"

template <typename real_t>
class eqs_todo_bulk_ode : public eqs_todo_bulk<real_t> 
{
  //// inheriting constructor (TODO: not suported by gcc as of 2012-05)
  //using eqs_todo_bulk<real_t>::eqs_todo_bulk;

  bool opt_cevp, opt_revp;
  public: eqs_todo_bulk_ode(const grd<real_t> &grid, map<enum eqs_todo_bulk<real_t>::processes, bool> opts);

  // RHS of the ODE to be solved for saturation adjustment 
  private: class rhs;

  private: void condevap(
    const mtx::arr<real_t> &rhod,
    mtx::arr<real_t> &rhod_th,
    mtx::arr<real_t> &rhod_rv,
    mtx::arr<real_t> &rhod_rl,
    mtx::arr<real_t> &rhod_rr,
    const quantity<si::time, real_t> dt
  );   
};
