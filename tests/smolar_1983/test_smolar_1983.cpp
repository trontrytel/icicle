/** @file
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 *  @section DESCRIPTION
 *  @section RESULTS
 *
 *    \image html "../../tests/smolar_1983/fig1.svg"
 */

#include <cstdlib>

int main() 
{ 
  int status = system("gdl -quiet -e test_smolar_1983");
  exit(status == 0 ? EXIT_SUCCESS : EXIT_FAILURE);
}
