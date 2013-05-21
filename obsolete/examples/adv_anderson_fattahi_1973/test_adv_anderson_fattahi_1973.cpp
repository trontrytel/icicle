/** @file
 *  @example adv_anderson_fattahi_1973/test_test_adv_anderson_fattahi_1973.cpp
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @section DESCRIPTION
      Rotating cone test from @copydetails Anderson_and_Fattahi_1974 J. Atmos. Sci.
 *  @section RESULTS
 *    \include "adv_anderson_fattahi_1973/test_adv_anderson_fattahi_1973.py"
 *    \image html "../../tests/adv_anderson_fattahi_1973/fig-leapfrog-2d.gif"
 *    \image html "../../tests/adv_anderson_fattahi_1973/fig-upstream-2d.gif"
 *    \image html "../../tests/adv_anderson_fattahi_1973/fig-mpdata-2d.gif"
 */

#include <cstdlib>

int main() 
{ 
  int status = system("python test_adv_anderson_fattahi_1973.py");
  exit(status == 0 ? EXIT_SUCCESS : EXIT_FAILURE);
}
