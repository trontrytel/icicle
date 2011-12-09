/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 *  @section DESCRIPTION
 *    A facility for indexing Blitz 3D arrays in a generic way.
 *    If a method has a template "idx" and then uses Array(idx(i,j,k))
 *    indexing, it can be called with idx=idx_ijk, idx=idx_jki
 *    or idx=idx_kij. Consequently if some calculations are done
 *    in the same way in all three directions, most of the code
 *    may be written once and just called thrice, each time wth
 *    a different value of the template.
 */
#ifndef IDX_HPP
#  define IDX_HPP

#  include "common.hpp" // Blitz includes

class idx_ijk : public RectDomain<3>
{
  public: idx_ijk(const Range &i, const Range &j, const Range &k) 
    : RectDomain<3>(TinyVector<Range, 3>(i,j,k)) 
  { } 
};

class idx_jki : public RectDomain<3>
{
  public: idx_jki(const Range &j, const Range &k, const Range &i) 
    : RectDomain<3>(TinyVector<Range, 3>(i,j,k)) 
  { } 
};

class idx_kij : public RectDomain<3>
{
  public: idx_kij(const Range &k, const Range &i, const Range &j) 
    : RectDomain<3>(TinyVector<Range, 3>(i,j,k)) 
  { } 
};

#endif
