/** @file
 *  @example adv_isolines/test_adv_isolines.cpp
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date February 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @section DESCRIPTION
      Truncation error test from Smolarkiewicz & Grabowski 1990 (J.Comp.Phys.,86,355-375)
 *  @section RESULTS
      Running icicle + calculating truncation error:
 *    \include "adv_isolines/isolines_1d.py"
      Plot:
 *    \include "adv_isolines/plot_iso.py"
      Results for MPDATA iord=2:
 *    \image html "../../tests/adv_isolines/fig-mpdata-iord2.svg"
 */

#include <cstdlib>

int main() 
{ 
  int status;

  status = system("python isolines_1d.py");
  exit(status == 0 ? EXIT_SUCCESS : EXIT_FAILURE);

  status = system("python plot_iso.py");
  exit(status == 0 ? EXIT_SUCCESS : EXIT_FAILURE);
}
