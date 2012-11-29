/** @file
 *  @example phc_kappa_koehler/phc_kappa_koehler_test.cpp
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date September 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @section DESCRIPTION
 *    Tests the implementation of the kappa-Koehler parametersiation
 *    following @copydetails Petters_and_Kreidenweis_2007 
 *    by checking if a drop growth rate would be equal to zero for 
 *    the equilibrium radius (i.e. if water activity would equal the RH)
 */

#include "../../src/phc/phc_kappa_koehler.hpp"
#include <list>
using std::list;

typedef double real_t;

int main()
{
  real_t eps = 1e-7; 

  // testing for several different kappas
  for (real_t &kappa : list<real_t>({2e-9, 2e-4, 2e-3, 2e-2, 2e-1, 2e0}))
  {
    // testing for a range of dry radii
    for (real_t &rd : list<real_t>({1e-9, 1e-8, 1e-7, 1e-6, 1e-5}))
    {
      // testing for a range of RH values
      for (real_t &vap_ratio : list<real_t>({.5, .6, .7, .8, .9, .999})
      )
      {
        quantity<si::dimensionless, real_t> aw = phc::kappa::a_w<real_t>(
          phc::kappa::rw3_eq_nokelvin<real_t>(pow<3>(rd * si::metres), kappa, vap_ratio),
          pow<3>(rd * si::metres), 
          kappa
        );
        if (std::abs(aw - vap_ratio) > eps) 
        {
          std::cerr << "aw - RH = " << aw - vap_ratio << std::endl;
          std::cerr << "for kappa = " << kappa << ", rd = " << rd << ", RH= " << vap_ratio << std::endl;
          exit(1);
        }
      }
    }
  }
}
