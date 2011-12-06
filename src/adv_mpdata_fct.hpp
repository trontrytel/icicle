/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 *  @section DESCRIPTION
 *    Flux Corrected Transport aka non-oscillatory aka monotonic
 *    option for MPDATA
 */
#ifndef ADV_MPDATA_FCT_HPP
#  define ADV_MPDATA_FCT_HPP

#  include "adv_mpdata.hpp"

template <typename real_t> 
class adv_mpdata_fct : public adv_mpdata<real_t> 
{
  public: adv_mpdata_fct(grd_arakawa_c_lorenz<real_t> *grid, int iord) 
    : adv_mpdata<real_t>(grid, iord, true)
  {
    // TODO: zaalokowac tmp1 i tmp2
  }

/*
  public: void op3D(
    Array<real_t, 3> *psi_ijk[],
    Array<real_t, 3> *psi_jki[],
    Array<real_t, 3> *psi_kij[],
    const Range &i,
    const Range &j,
    const Range &k,
    const int n, const int step,
    arr<real_t> &Cx,
    arr<real_t> &Cy,
    arr<real_t> &Cz
  )
  {
    adv_mpdata<real_t>::op3D(psi_ijk, psi_jki, psi_kij, i, j, k, n, step, Cx, Cy, Cz);
    if (step > 1)
    {
      // calculating Cx_mon, Cy_mon, Cz_mon
//      adv_mpdata<real_t>::op(psi_ijk, psi_jki, psi_kij, i, j, k, n, step, Cx_mon, Cy_mon, Cz_mon);
    }
  }
*/

};
#endif
