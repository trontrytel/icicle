/** @file
 *  @example eqs_harmonic_oscillator/test_eqs_harmonic_oscillator.cpp
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date March 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @section DESCRIPTION
      test for simplecticity from Smolarkiewicz 2006 Int. J. Numer. Meth. Fluids 
 *  @section RESULTS
 *    \image html "../../tests/eqs_harmonic_oscillator/harmonic_osc.gif"
 *    \include "eqs_harmonic_oscillator/test_harmonic_oscillator.py"
 */

#include <cstdlib>
extern "C" {
#include <unistd.h>
}

int main() 
{ 
  unlink("out.nc");
  int status = system("python test_harmonic_oscillator.py");
  exit(status == 0 ? EXIT_SUCCESS : EXIT_FAILURE);
}
