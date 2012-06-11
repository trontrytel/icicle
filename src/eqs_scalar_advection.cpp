/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date January 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @brief contains definition of the eqs_scalar_advection class - the simplest equation system
 */

#include "cfg.hpp" 
#include "eqs_scalar_advection.hpp"

/// @brief the simplest equation system consisting of a single homogeneous transport equation
template <typename real_t>
eqs_scalar_advection<real_t>::eqs_scalar_advection()
{
  this->sys.push_back(new struct eqs<real_t>::gte({
    "psi", "the transported scalar field", 
    this->quan2str(quantity<si::dimensionless, real_t>(1.))
  }));
}

// explicit instantiations
#if defined(USE_FLOAT)
template class eqs_scalar_advection<float>;
#endif
#if defined(USE_DOUBLE)
template class eqs_scalar_advection<double>;
#endif
#if defined(USE_LDOUBLE)
template class eqs_scalar_advection<long double>;
#endif