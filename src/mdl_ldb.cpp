/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */

#include "cfg.hpp"
#include "cmn.hpp"
#include "mdl.hpp"

void mdl_ldb(const po::variables_map &vm, const string &options) 
{
  mdl<long double>(vm, options);
}
