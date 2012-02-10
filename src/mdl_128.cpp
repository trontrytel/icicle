/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date February 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#include "cfg.hpp"
#include "cmn.hpp"
#include "mdl.hpp"

void mdl_128(const po::variables_map &vm, const string &options) 
{
  // TODO: error if not available with the current compiler...
  mdl<__float128>(vm, options);
}
