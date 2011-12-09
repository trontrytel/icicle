/** @file
 *  @example smolar_1983/test_smolar_1983.cpp
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 *  @section DESCRIPTION
 *  @section RESULTS
 *    \include "smolar_1983/test_smolar_1983.py"
 *    \image html "../../tests/smolar_1983/fig-leapfrog.gif"
 *    \image html "../../tests/smolar_1983/fig-upstream.gif"
 *    \image html "../../tests/smolar_1983/fig-mpdata.gif"
 */

#include <cstdlib>

int main() 
{ 
  int status = system("python test_smolar_1983.py");
  exit(status == 0 ? EXIT_SUCCESS : EXIT_FAILURE);
}
