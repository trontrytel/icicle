/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date February 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#include "cfg.hpp"
#include "vel_momeq_extrapol.hpp" 
#include "grd_carthesian.hpp" 

template <typename real_t>
void vel_momeq_extrapol<real_t>::populate_courant_fields(int nm0, int nm1,
    mtx::arr<real_t> *Cx, 
    mtx::arr<real_t> *Cy, 
    mtx::arr<real_t> *Cz, 
    quantity<si::time, real_t> dt,
    mtx::arr<real_t> *Qx[],
    mtx::arr<real_t> *Qy[],
    mtx::arr<real_t> *Qz[]
) const
{
    {
      mtx::arr<real_t> *C[] = { Cx, Cy, Cz };
      mtx::arr<real_t> **Q[] = { Qx, Qy, Qz };
      for (int xyz = 0; xyz < 3; ++xyz)
      {
        if (Q[xyz] != NULL) 
        {
          const mtx::idx 
            imph_j_k(mtx::idx_ijk(C[xyz]->i - grd_carthesian<real_t>::p_half /* sic! */, C[xyz]->j, C[xyz]->k)),
            ipmh_j_k(mtx::idx_ijk(C[xyz]->i + grd_carthesian<real_t>::m_half /* sic! */, C[xyz]->j, C[xyz]->k));

          (*C[xyz])(C[xyz]->ijk) = extrapol(
            (*Q[xyz][nm0])(imph_j_k),
            (*Q[xyz][nm0])(ipmh_j_k),
            (*Q[xyz][nm1])(imph_j_k),
            (*Q[xyz][nm1])(ipmh_j_k)
          );
        }
      }
    }
    if (Qx != NULL) (*Cx)(Cx->ijk) *= (dt/si::seconds) / (grid.dx()/si::metres);
    if (Qy != NULL) (*Cy)(Cy->ijk) *= (dt/si::seconds) / (grid.dy()/si::metres);
    if (Qz != NULL) (*Cz)(Cz->ijk) *= (dt/si::seconds) / (grid.dz()/si::metres);
}

// explicit instantiations
#if defined(USE_FLOAT)
template class vel_momeq_extrapol<float>;
#endif
#if defined(USE_DOUBLE)
template class vel_momeq_extrapol<double>;
#endif
#if defined(USE_LDOUBLE)
template class vel_momeq_extrapol<long double>;
#endif
