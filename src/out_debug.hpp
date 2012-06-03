/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef OUT_DEBUG_HPP
#  define OUT_DEBUG_HPP

#  include "out.hpp"

extern "C" {
#  include <unistd.h>
}

template <typename real_t>
class out_debug : public out<real_t>
{
  public: virtual void record(
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
      << endl
      << psi(ijk) 
      << endl;
    std::cerr << tmp.str(); // non-buffered?
  }
};

#endif
