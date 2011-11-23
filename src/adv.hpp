/** @file
 *  @author Sylwester Arabas (University of Warsaw) <slayoo@igf.fuw.edu.pl>
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef ADV_HPP
#  define ADV_HPP

#  include "config.hpp" // USE_* defines
#  include "common.hpp" // root class, error reporting, ...
#  include "grd.hpp" // m_half, p_half, ...

template <class unit, typename real_t>
class adv : root
{
  public: virtual const int stencil_extent() = 0;
  public: virtual const int time_levels() = 0;
  public: virtual const int num_steps() = 0;

  public: virtual void op_1D(Array<quantity<unit, real_t>, 3> *psi[], 
    const Range &i, 
    const Range &j, 
    const Range &k, 
    const int n, const int step,
    const Array<quantity<si::dimensionless, real_t>, 3> &Cx,
    const Array<quantity<si::dimensionless, real_t>, 3> &Cy,
    const Array<quantity<si::dimensionless, real_t>, 3> &Cz
  ) = 0;
};
#endif
