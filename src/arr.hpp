/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef ARR_HPP
#  define ARR_HPP

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

/**   @class idx
 *    A facility for indexing Blitz 3D arrays in a generic way.
 *    If a method has a template "idx" and then uses arr(idx(i,j,k))
 *    indexing, it can be called with idx=idx_ijk, idx=idx_jki
 *    or idx=idx_kij. Consequently if some calculations are done
 *    in the same way in all three directions, most of the code
 *    may be written once and just called thrice, each time wth
 *    a different value of the template.
 */
class idx : public RectDomain<3> 
{
  public: idx(const TinyVector<Range, 3> &tv)
    : RectDomain<3>(tv)
  { }
};

class idx_ijk : public idx
{
  public: idx_ijk(const Range &i, const Range &j, const Range &k) 
    : idx(TinyVector<Range, 3>(i,j,k)) 
  { } 
};

class idx_jki : public idx
{
  public: idx_jki(const Range &j, const Range &k, const Range &i) 
    : idx(TinyVector<Range, 3>(i,j,k)) 
  { } 
};

class idx_kij : public idx
{
  public: idx_kij(const Range &k, const Range &i, const Range &j) 
    : idx(TinyVector<Range, 3>(i,j,k)) 
  { } 
};

#endif
