/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef OUT_HPP
#  define OUT_HPP

#  include "cmn.hpp" // root class
#  include "mtx.hpp" // Blitz includes

template <typename real_t>
class out : root
{
  public: virtual void record( // TODO: t should be in seconds
    mtx::arr<real_t> *psi, 
    const mtx::rng &i, const mtx::rng &j, const mtx::rng &k, const unsigned long t
  ) = 0;
};

#endif
