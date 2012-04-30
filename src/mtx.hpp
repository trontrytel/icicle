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
using blitz::where;

#  define mtx_expr_1arg_macro(name,a1,expr) \
  template<class t1> \
  static auto name(const t1 &a1) \
  -> decltype(blitz::safeToReturn(expr)) \
  { \
    assert(isfinite(sum(a1))); \
    return blitz::safeToReturn(expr); \
  } 
#  define mtx_expr_2arg_macro(name,a1,a2,expr) \
  template<class t1, class t2> \
  static auto name(const t1 &a1, const t2 &a2) \
  -> decltype(blitz::safeToReturn(expr)) \
  { \
    assert(isfinite(sum(a1))); \
    assert(isfinite(sum(a2))); \
    return blitz::safeToReturn(expr); \
  } 
#  define mtx_expr_3arg_macro(name, a1, a2, a3, expr) \
  template<class t1, class t2, class t3> \
  static auto name(const t1 &a1, const t2 &a2, const t3 &a3) \
  -> decltype(blitz::safeToReturn(expr)) \
  { \
    assert(isfinite(sum(a1))); \
    assert(isfinite(sum(a2))); \
    assert(isfinite(sum(a3))); \
    return blitz::safeToReturn(expr); \
  }
#  define mtx_expr_4arg_macro(name, a1, a2, a3, a4, expr) \
  template<class t1, class t2, class t3, class t4> \
  static auto name(const t1 &a1, const t2 &a2, const t3 &a3, const t4 &a4) \
  -> decltype(blitz::safeToReturn(expr)) \
  { \
    assert(isfinite(sum(a1))); \
    assert(isfinite(sum(a2))); \
    assert(isfinite(sum(a3))); \
    assert(isfinite(sum(a4))); \
    return blitz::safeToReturn(expr); \
  }
#  define mtx_expr_5arg_macro(name, a1, a2, a3, a4, a5, expr) \
  template<class t1, class t2, class t3, class t4, class t5> \
  static auto name(const t1 &a1, const t2 &a2, const t3 &a3, const t4 &a4, const t5 &a5) \
  -> decltype(blitz::safeToReturn(expr)) \
  { \
    assert(isfinite(sum(a1))); \
    assert(isfinite(sum(a2))); \
    assert(isfinite(sum(a3))); \
    assert(isfinite(sum(a4))); \
    assert(isfinite(sum(a5))); \
    return blitz::safeToReturn(expr); \
  }
#  define mtx_expr_6arg_macro(name, a1, a2, a3, a4, a5, a6, expr) \
  template<class t1, class t2, class t3, class t4, class t5, class t6> \
  static auto name(const t1 &a1, const t2 &a2, const t3 &a3, const t4 &a4, const t5 &a5, const t6 &a6) \
  -> decltype(blitz::safeToReturn(expr)) \
  { \
    assert(isfinite(sum(a1))); \
    assert(isfinite(sum(a2))); \
    assert(isfinite(sum(a3))); \
    assert(isfinite(sum(a4))); \
    assert(isfinite(sum(a5))); \
    assert(isfinite(sum(a6))); \
    return blitz::safeToReturn(expr); \
  }
#  define mtx_expr_11arg_macro(name, a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11, expr) \
  template<class t1, class t2, class t3, class t4, class t5, class t6, class t7, class t8, class t9, class t10, class t11> \
  static auto name(const t1 &a1, const t2 &a2, const t3 &a3, const t4 &a4, const t5 &a5, const t6 &a6, const t7 &a7, const t8 &a8, const t9 &a9, const t10 &a10, const t11 &a11) \
  -> decltype(blitz::safeToReturn(expr)) \
  { \
    assert(isfinite(sum(a1))); \
    assert(isfinite(sum(a2))); \
    assert(isfinite(sum(a3))); \
    assert(isfinite(sum(a4))); \
    assert(isfinite(sum(a5))); \
    assert(isfinite(sum(a6))); \
    assert(isfinite(sum(a7))); \
    assert(isfinite(sum(a8))); \
    assert(isfinite(sum(a9))); \
    assert(isfinite(sum(a10))); \
    assert(isfinite(sum(a11))); \
    return blitz::safeToReturn(expr); \
  }

#  undef guard

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

///   @brief A facility for indexing Blitz 3D blitz::Arrays in a generic way.
/**
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

  /// @brief non-permuted version of idx
  struct idx_ijk : idx
  {
    idx_ijk(const blitz::Range &i, const blitz::Range &j, const blitz::Range &k) 
      : idx(blitz::TinyVector<blitz::Range, 3>(i,j,k)) 
    { } 
    idx_ijk(const blitz::Range &i, const blitz::Range &j, const int k = 0) 
      : idx(blitz::TinyVector<blitz::Range, 3>(i,j,blitz::Range(k,k))) 
    { } 
    idx_ijk(const blitz::Range &i, const int j = 0, const int k = 0) 
      : idx(blitz::TinyVector<blitz::Range, 3>(
        i,
        blitz::Range(j,j),
        blitz::Range(k,k)
      )) 
    { } 
    idx_ijk(const int i, const int j = 0, const int k = 0) 
      : idx(blitz::TinyVector<blitz::Range, 3>(
        blitz::Range(i,i),
        blitz::Range(j,j),
        blitz::Range(k,k)
      )) 
    { } 
  };

  /// @brief ijk->jki permuted version of idx
  struct idx_jki : idx
  {
    idx_jki(const blitz::Range &j, const blitz::Range &k, const blitz::Range &i) 
      : idx(blitz::TinyVector<blitz::Range, 3>(i,j,k)) 
    { } 
  };

  /// @brief ijk->kij permuted version of idx
  struct idx_kij : idx
  {
    idx_kij(const blitz::Range &k, const blitz::Range &i, const blitz::Range &j) 
      : idx(blitz::TinyVector<blitz::Range, 3>(i,j,k)) 
    { } 
  };

  /// @brief 3D array representation - a thin wrapper over the Array class from the Blitz++ library
  template <typename real_t>
  struct arr : blitz::Array<real_t, 3>
  {
    typedef blitz::Array<real_t, 3> type;

    idx ijk;
    rng i, j, k;

    void fill_with_nans() 
    { 
      (*this)(ijk) = nan<real_t>(); 
    }

    // ctor
    arr(const idx &ijk) :
      blitz::Array<real_t, 3>(
        blitz::Range(ijk.lbound(0), ijk.ubound(0)),
        blitz::Range(ijk.lbound(1), ijk.ubound(1)),
        blitz::Range(ijk.lbound(2), ijk.ubound(2)),
        blitz::ColumnMajorArray<3>()
      ),
      i(blitz::Range(ijk.lbound(0), ijk.ubound(0))),
      j(blitz::Range(ijk.lbound(1), ijk.ubound(1))),
      k(blitz::Range(ijk.lbound(2), ijk.ubound(2))),
      ijk(ijk)
    {
#  ifndef NDEBUG
      fill_with_nans();
#  endif
    }
  };
};
#endif
