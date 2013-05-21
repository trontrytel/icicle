/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date April-June 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @brief a set of enums used in sdm
 */
#pragma once

namespace sdm 
{
  enum thrust_device_systems {openmp, cuda};
  enum chem_gas {gSO2, gO3, gH2O2};
  enum chem_aq {H, OH, SO2, O3, H2O2, HSO3, SO3, S_VI ,HSO4, SO4};
  enum processes {adve, cond, sedi, coal, chem};
  enum ode_algos {euler, mmid, rk4};
  enum xi_dfntns {id, ln, p2, p3};
};
