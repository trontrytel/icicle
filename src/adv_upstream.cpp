/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @brief definition of the adv_upstream class with an implementation of the upstream advection scheme
 */
#include "cfg.hpp"
#include "adv_upstream.hpp"

template <typename real_t>
typename adv<real_t>::op3D *adv_upstream<real_t>::factory(
  const mtx::idx &ijk,
  mtx::arr<real_t> **, 
  mtx::arr<real_t> **, 
  bool
) const
{
  return new op3D(ijk);
}

// explicit instantiations
#if defined(USE_FLOAT)
template class adv_upstream<float>;
#endif
#if defined(USE_DOUBLE)
template class adv_upstream<double>;
#endif
#if defined(USE_LDOUBLE)
template class adv_upstream<long double>;
#endif
