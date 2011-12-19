/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef ARR_HPP
#  define ARR_HPP

#  include "config.hpp"
#  if defined(USE_BOOST_THREAD) || defined(_OPENMP)
#    define BZ_THREADSAFE
#  endif
#  include <blitz/array.h>
#  include <blitz/numinquire.h>
using namespace blitz;

template <typename real_t>
class arr : public Array<real_t, 3>
{
  public: arr(const Range &i, const Range &j, const Range &k) 
    : Array<real_t, 3>(i, j, k)
  { 
    (*this)(i,j,k) = has_signalling_NaN(real_t(0)) 
      ? signalling_NaN(real_t(0)) 
      : quiet_NaN(real_t(0));
  } 
};

#endif
