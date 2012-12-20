/** @file
 *  @example eqs_shallow_water/test_eqs_shallow_water.cpp
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date March 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @section DESCRIPTION
 *    \image html "../../tests/eqs_shallow_water/fig1.svg"
 *    \include "eqs_shallow_water/test_eqs_shallow_water.py"
 */

#include <cstdlib>
extern "C" {
#include <unistd.h>
}

int main() 
{ 
  unlink("out.nc");
  int status = system("python test_eqs_shallow_water.py");
  exit(status == 0 ? EXIT_SUCCESS : EXIT_FAILURE);
}
