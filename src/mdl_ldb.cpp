/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */

#include "config.hpp"
#include "common.hpp"
#include "mdl.hpp"

void mdl_ldb(const po::variables_map& vm) 
{
  mdl<long double>(vm);
}
