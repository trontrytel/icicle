/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef OUT_HPP
#  define OUT_HPP

#  include "config.hpp" // USE_* defines
#  include "common.hpp" // root class, error reporting

template <typename real_t>
class out : root
{
  public: virtual void record( // TODO: t should be in seconds
    Array<real_t, 3> **psi, const int n, 
    const Range &i, const Range &j, const Range &k, const unsigned long t
  ) = 0;
};

#endif
