/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date February 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#pragma once
#include "vel_momeq.hpp" 
#include "grd.hpp"

template <typename real_t>
class vel_momeq_extrapol : public vel_momeq<real_t>
{
  private: const grd<real_t> &grid;

  public: vel_momeq_extrapol(const grd<real_t> &grid)
    : grid(grid)
  { }

  // eq. 34b in Smolarkiewicz & Margolin 1998
  // .25 = .5 * .5 (stems from i +/- 1/2 averaging)
  private: mtx_expr_4arg_macro(extrapol, nm0_il, nm0_ir, nm1_il, nm1_ir,
    real_t(.25) * (real_t(3) * (nm0_il + nm0_ir) - (nm1_il + nm1_ir)) 
  )

  public: void populate_courant_fields(int nm0, int nm1,
    mtx::arr<real_t> *Cx, 
    mtx::arr<real_t> *Cy, 
    mtx::arr<real_t> *Cz, 
    quantity<si::time, real_t> dt,
    mtx::arr<real_t> *Qx[],
    mtx::arr<real_t> *Qy[],
    mtx::arr<real_t> *Qz[]
  ) const;
};
