/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011 - March 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @brief definition of the adv_mpdata class with an implementation of the MPDATA advection scheme
 */

#include "adv_mpdata.hpp"

template <typename real_t>
adv_mpdata<real_t>::adv_mpdata(int iord, bool cross_terms, bool third_order) // TODO: enums?
  : iord(iord), cross_terms(cross_terms), third_order(third_order), adv_upstream<real_t>()
{
  if (iord <= 0) error_macro("iord (the number of iterations) must be > 0")
  if (iord < 3 && third_order) warning_macro("third-order accuracy needs iord >= 3")
}

template <typename real_t>
typename adv<real_t>::op3D *adv_mpdata<real_t>::factory(
  const mtx::idx &ijk,
  mtx::arr<real_t> **tmp_s,
  mtx::arr<real_t> **tmp_v,
  bool positive_definite
) const
{
  if (positive_definite)
    return new op3D<aon_nil>(ijk, tmp_s, tmp_v, cross_terms, iord, third_order);
  else
    return new op3D<aon_abs>(ijk, tmp_s, tmp_v, cross_terms, iord, third_order);
}

// explicit instantiations
#define ICICLE_INSTANTIATE_CLASS adv_mpdata
#include "cmn/cmn_instant.hpp"
