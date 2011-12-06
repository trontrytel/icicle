/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
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
