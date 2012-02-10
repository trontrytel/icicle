/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef MTX_HPP
#  define MTX_HPP

#  if defined(USE_BOOST_THREAD) || defined(_OPENMP)
#    define BZ_THREADSAFE
#  endif
#  include <blitz/array.h>
#  include <blitz/numinquire.h>

namespace mtx
{
  const int i = blitz::firstRank;
  const int j = blitz::secondRank;
  const int k = blitz::thirdRank;

  typedef blitz::Range rng;

  template <typename real_t>
  real_t nan()
  {
    return blitz::has_signalling_NaN(real_t(0)) 
      ? blitz::signalling_NaN(real_t(0)) 
      : blitz::quiet_NaN(real_t(0));
  }

  template <typename real_t>
  real_t eps()
  {
    return blitz::epsilon(real_t(0));
  }

  template <typename arr> void cycle(arr &a, arr &b) 
  { 
    blitz::cycleArrays(a, b); 
  }

  template <typename arr> void cycle(arr &a, arr &b, arr &c) 
  { 
    blitz::cycleArrays(a, b, c); 
  }

/**   @class idx
 *    A facility for indexing Blitz 3D blitz::Arrays in a generic way.
 *    If a method has a template "idx" and then uses arr(idx(i,j,k))
 *    indexing, it can be called with idx=idx_ijk, idx=idx_jki
 *    or idx=idx_kij. Consequently if some calculations are done
 *    in the same way in all three directions, most of the code
 *    may be written once and just called thrice, each time wth
 *    a different value of the template.
 */
  struct idx : blitz::RectDomain<3> 
  {
    rng i, j, k;
    bool i_spans, j_spans, k_spans;

    idx(const blitz::TinyVector<blitz::Range, 3> &tv)
      : blitz::RectDomain<3>(tv), 
        i(tv[mtx::i]), j(tv[mtx::j]), k(tv[mtx::k]), 
        i_spans(i.first() != i.last()),
        j_spans(j.first() != j.last()),
        k_spans(k.first() != k.last())
    { }
  };

  struct idx_ijk : idx
  {
    idx_ijk(const blitz::Range &i, const blitz::Range &j, const blitz::Range &k) 
      : idx(blitz::TinyVector<blitz::Range, 3>(i,j,k)) 
    { } 
  };

  struct idx_jki : idx
  {
    idx_jki(const blitz::Range &j, const blitz::Range &k, const blitz::Range &i) 
      : idx(blitz::TinyVector<blitz::Range, 3>(i,j,k)) 
    { } 
  };

  struct idx_kij : idx
  {
    idx_kij(const blitz::Range &k, const blitz::Range &i, const blitz::Range &j) 
      : idx(blitz::TinyVector<blitz::Range, 3>(i,j,k)) 
    { } 
  };

  template <typename real_t>
  struct arr : blitz::Array<real_t, 3>
  {
    typedef blitz::Array<real_t, 3> type;

    arr(const idx &ijk)
      : blitz::Array<real_t, 3>(
        blitz::Range(ijk.lbound(0), ijk.ubound(0)), 
        blitz::Range(ijk.lbound(1), ijk.ubound(1)), 
        blitz::Range(ijk.lbound(2), ijk.ubound(2))
      )
    {
#  ifndef NDEBUG
      (*this)(ijk) = nan<real_t>();
#  endif
    }
  };
};

#endif
