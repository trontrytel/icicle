/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef OUT_DEBUG_HPP
#  define OUT_DEBUG_HPP

#  include "out.hpp"

template <typename real_t>
class out_debug : public out<real_t>
{
  public: virtual void record(
    Array<real_t, 3> **psi, const int n, 
    const Range &i, const Range &j, const Range &k, const unsigned long t
  ) 
  {
    cerr // non-buffered?
      << "[" << i << "," << j << "," << k << "] @ t/dt=" << t
      << endl
      << (*psi[n])(i, j, k) 
      << endl;
  }
};

#endif
