/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#include "cfg.hpp"
#include "out_debug.hpp"

extern "C" {
#  include <unistd.h>
}

template <typename real_t> 
void out_debug<real_t>::record(
    const string &name,
    const mtx::arr<real_t> &psi,
    const mtx::idx &ijk, 
    const unsigned long t
)
{
    assert(name == "psi"); // TODO
    std::ostringstream tmp;
    tmp
      << "{pid: " << getpid() << "} :" 
      << "[" 
        << ijk.lbound(0) << "..." << ijk.ubound(0) << "," 
        << ijk.lbound(1) << "..." << ijk.ubound(1) << "," 
        << ijk.lbound(2) << "..." << ijk.ubound(2) << "," 
      << "] @ t/dt=" << t
      << std::endl
      << psi(ijk) 
      << std::endl;
    std::cerr << tmp.str(); // non-buffered?
}

// explicit instantiations
#if defined(USE_FLOAT)
template class out_debug<float>;
#endif
#if defined(USE_DOUBLE)
template class out_debug<double>;
#endif
#if defined(USE_LDOUBLE)
template class out_debug<long double>;
#endif
