/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date December 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 *  @section DESCRIPTION
 *    A simple container for storing simulation set-up elements, i.e.
 *    the advection scheme, the output type, the velocity field etc
 */
#ifndef STP_HPP
#  define SPT_HPP

#  include "adv.hpp"
#  include "out.hpp"
#  include "vel.hpp"
#  include "arr.hpp"
#  include "grd.hpp"
#  include "ini.hpp"

template <typename real_t>
class stp : root
{
  // that should probably be the only place where public fields are used
  public: adv<real_t> *fllbck, *advsch;
  public: out<real_t> *output;
  public: vel<real_t> *velocity;
  public: ini<real_t> *intcond;
  public: grd<real_t> *grid;

  public: stp(
    adv<real_t> *fllbck, adv<real_t> *advsch, 
    out<real_t> *output, 
    vel<real_t> *velocity,
    ini<real_t> *intcond,
    grd<real_t> *grid
  ) 
    : fllbck(fllbck), advsch(advsch), 
      output(output), 
      velocity(velocity), 
      intcond(intcond), 
      grid(grid)
  {}
};

#endif
