/** @file
 *  @example eqs_isentropic_halny/test_eqs_isentropic_halny.cpp
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date March 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @section DESCRIPTION
 *    \image html "../../tests/eqs_isentropic_halny/halny.gif"
 *    \include "eqs_isentropic_halny/test_eqs_isentropic_halny.py"
 */

#include <cstdlib>
extern "C" {
#include <unistd.h>
}

int main() 
{ 
  unlink("out.nc");
  int status = system("python test_eqs_isentropic_halny.py");
  exit(status == 0 ? EXIT_SUCCESS : EXIT_FAILURE);
}
