/** @file
 *  @example mayavi/test_mayavi.cpp
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @section DESCRIPTION
 *  @section RESULTS
 *    \include "mayavi/test_mayavi.py"
 *    \image html "../../tests/mayavi/fig-leapfrog-3d.gif"
 *    \image html "../../tests/mayavi/fig-upstream-3d.gif"
 *    \image html "../../tests/mayavi/fig-mpdata-3d.gif"
 */

#include <cstdlib>

int main() 
{ 
  int status = system("python test_mayavi.py");
  exit(status == 0 ? EXIT_SUCCESS : EXIT_FAILURE);
}
