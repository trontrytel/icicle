/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date September 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @brief explicit template instantiation automation
 */

#pragma once

#include "../cfg/cfg_types.hpp"

#if !defined(ICICLE_INSTANTIATE_CLASS)
#  error "ICICLE_INSTANTIATE_CLASS not defined"
#endif

#if !defined(ICICLE_INSTANTIATE_PARAM)

#  if defined(USE_FLOAT)
template class ICICLE_INSTANTIATE_CLASS<float>;
#  endif
#  if defined(USE_DOUBLE)
template class ICICLE_INSTANTIATE_CLASS<double>;
#  endif
#  if defined(USE_LDOUBLE)
template class ICICLE_INSTANTIATE_CLASS<long double>;
#  endif
#  if defined(USE_FLOAT128)
template class ICICLE_INSTANTIATE_CLASS<__float128>;
#  endif

#else

#  if defined(ICICLE_INSTANTIATE_PARAM_PARAM)

#    if defined(USE_FLOAT)
template class ICICLE_INSTANTIATE_CLASS<float, ICICLE_INSTANTIATE_PARAM<float>>;
#    endif
#    if defined(USE_DOUBLE)
template class ICICLE_INSTANTIATE_CLASS<double, ICICLE_INSTANTIATE_PARAM<double>>;
#    endif
#    if defined(USE_LDOUBLE)
template class ICICLE_INSTANTIATE_CLASS<long double, ICICLE_INSTANTIATE_PARAM<long double>>;
#    endif
#    if defined(USE_FLOAT128)
template class ICICLE_INSTANTIATE_CLASS<__float128, ICICLE_INSTANTIATE_PARAM<__float128>>;
#    endif

#  else

#    if defined(USE_FLOAT)
template class ICICLE_INSTANTIATE_CLASS<float, ICICLE_INSTANTIATE_PARAM>;
#    endif
#    if defined(USE_DOUBLE)
template class ICICLE_INSTANTIATE_CLASS<double, ICICLE_INSTANTIATE_PARAM>;
#    endif
#    if defined(USE_LDOUBLE)
template class ICICLE_INSTANTIATE_CLASS<long double, ICICLE_INSTANTIATE_PARAM>;
#    endif
#    if defined(USE_FLOAT128)
template class ICICLE_INSTANTIATE_CLASS<__float128, ICICLE_INSTANTIATE_PARAM>;
#    endif

#  endif

#endif
