/** @file
 *  @author Sylwester Arabas (University of Warsaw) <slayoo@igf.fuw.edu.pl>
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef GRD_HPP
#  define GRD_HPP

#  include "config.hpp" // USE_* defines
#  include "common.hpp" // root class, error reporting

class grd : root
{
  public: static const int m_half = 0;
  public: static const int p_half = 1;
};
#endif
