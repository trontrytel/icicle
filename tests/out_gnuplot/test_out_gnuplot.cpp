/** @file
 *  @example gnuplot/test_gnuplot.cpp
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date December 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @section DESCRIPTION
 *    The out_gnuplot module (chosen with the --out=gnuplot option) 
 *      is designed to ease visualisation of 1D simulations in gnuplot.
 *    An example gnuplot script and the resultant SVG file are shown below:
 *    \include "out_gnuplot/test_out_gnuplot.gpi"
 *    \image html "../../tests/out_gnuplot/fig1.svg"
 */

#include <cstdlib>

int main() 
{ 
  if (0 != system("gnuplot test_out_gnuplot.gpi")) 
    exit(EXIT_FAILURE);
}
