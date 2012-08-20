/** @file
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date August 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 @  @brief contains definitions of eqs_todo_mm class - 2 moment parametrisation of warm rain microphysics
 */

#pragma once
#include "eqs_todo.hpp"

template <typename real_t>
class eqs_todo_mm : public eqs_todo<real_t> 
{
  // nested class (well... struct)
  protected: struct params : eqs_todo<real_t>::params
  {
    int idx_rhod_rl, idx_rhod_rr, idx_rhod_nl, idx_rhod_nr; // advected variables indices
  };

  protected: params par;

  // a container for storing options (i.e. which processes ar on/off)
  public: enum processes {act, cond, acc, autoc, self, turb, sedi};
  private: map<enum processes, bool> opts;

  // nested class
  private: template <bool fill_tmp> class activation;
  private: class cond_evap;


  // ctor
  public: eqs_todo_mm(
    const grd<real_t> &grid, 
    map<enum processes, bool> opts, 
    real_t mean_rd, 
    real_t sdev_rd,
    real_t n_tot,
    real_t chem_b
  );
};
