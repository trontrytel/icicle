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
//using namespace blitz; // TODO: to be removed!

typedef blitz::Range rng;

template <typename real_t>
class arr : public blitz::Array<real_t, 3>
{
  public: typedef blitz::Array<real_t, 3> arr_ret;

  public: arr(const blitz::Range &i, const blitz::Range &j, const blitz::Range &k) 
    : blitz::Array<real_t, 3>(i, j, k)
  { 
// TODO: move to solver, and make arr a typedef..
// TODO: arr::nan arr::eps ...?
    (*this)(i,j,k) = blitz::has_signalling_NaN(real_t(0)) 
      ? blitz::signalling_NaN(real_t(0)) 
      : blitz::quiet_NaN(real_t(0));
  } 
};

/**   @class idx
 *    A facility for indexing Blitz 3D blitz::Arrays in a generic way.
 *    If a method has a template "idx" and then uses arr(idx(i,j,k))
 *    indexing, it can be called with idx=idx_ijk, idx=idx_jki
 *    or idx=idx_kij. Consequently if some calculations are done
 *    in the same way in all three directions, most of the code
 *    may be written once and just called thrice, each time wth
 *    a different value of the template.
 */
class idx : public blitz::RectDomain<3> 
{
  public: idx(const blitz::TinyVector<blitz::Range, 3> &tv)
    : blitz::RectDomain<3>(tv)
  { }
};

class idx_ijk : public idx
{
  public: idx_ijk(const blitz::Range &i, const blitz::Range &j, const blitz::Range &k) 
    : idx(blitz::TinyVector<blitz::Range, 3>(i,j,k)) 
  { } 
};

class idx_jki : public idx
{
  public: idx_jki(const blitz::Range &j, const blitz::Range &k, const blitz::Range &i) 
    : idx(blitz::TinyVector<blitz::Range, 3>(i,j,k)) 
  { } 
};

class idx_kij : public idx
{
  public: idx_kij(const blitz::Range &k, const blitz::Range &i, const blitz::Range &j) 
    : idx(blitz::TinyVector<blitz::Range, 3>(i,j,k)) 
  { } 
};

#endif
