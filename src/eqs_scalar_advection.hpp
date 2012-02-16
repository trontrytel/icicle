/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date January 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef EQS_SCALAR_ADVECTION_HPP
#  define EQS_SCALAR_ADVECTION_HPP

#  include "cmn.hpp" // root class, error reporting
#  include "eqs.hpp"

template <typename real_t>
class eqs_scalar_advection : public eqs<real_t>
{
  private: vector<struct eqs<real_t>::gte> sys;
  public: vector<struct eqs<real_t>::gte> & system() { return sys; }

  public: eqs_scalar_advection()
  {
    {
      struct eqs<real_t>::gte e = {
        "psi", "the transported scalar field", 
        this->quan2str(quantity<si::dimensionless, real_t>(1.)),
        0,0,0 // TODO: some nicer syntax ...
      };
      sys.push_back(e);
    }
  }
};
#endif
