/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date February 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef VEL_MOMEQ_EXTRAPOL_HPP
#  define VEL_MOMEQ_EXTRAPOL_HPP

#  include "vel_momeq.hpp" 
#  include "grd_arakawa-c-lorenz.hpp"

template <typename real_t>
class vel_momeq_extrapol : public vel_momeq<real_t>
{
  private: grd_arakawa_c_lorenz<real_t> *grid;

  public: vel_momeq_extrapol(grd<real_t> *grid)
    : grid(dynamic_cast<grd_arakawa_c_lorenz<real_t>*>(grid))
  {
    if (grid == NULL) error_macro("vel_momeq_extrapol supports the Arakawa-C grid only")
  }

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
  )   
  {
    {
      mtx::arr<real_t> *C[] = { Cx, Cy, Cz };
      mtx::arr<real_t> **Q[] = { Qx, Qy, Qz };
      for (int xyz = 0; xyz < 3; ++xyz)
      {
        if (Q[xyz] != NULL) 
        {
          const mtx::idx 
            imph_j_k(mtx::idx_ijk(C[xyz]->i - grid->p_half /* sic! */, C[xyz]->j, C[xyz]->k)),
            ipmh_j_k(mtx::idx_ijk(C[xyz]->i + grid->m_half /* sic! */, C[xyz]->j, C[xyz]->k));

          (*C[xyz])(C[xyz]->ijk) = extrapol(
            (*Q[xyz][nm0])(imph_j_k),
            (*Q[xyz][nm0])(ipmh_j_k),
            (*Q[xyz][nm1])(imph_j_k),
            (*Q[xyz][nm1])(ipmh_j_k)
          );
        }
      }
    }
    if (Qx != NULL) (*Cx)(Cx->ijk) *= (dt/si::seconds) / (grid->dx()/si::metres);
    if (Qy != NULL) (*Cy)(Cy->ijk) *= (dt/si::seconds) / (grid->dy()/si::metres);
    if (Qz != NULL) (*Cz)(Cz->ijk) *= (dt/si::seconds) / (grid->dz()/si::metres);
  }
};
#endif
