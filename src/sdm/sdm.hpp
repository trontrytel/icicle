/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date April-June 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#pragma once

#include "../cmn/cmn_error.hpp"
#include "../phc/phc.hpp"
#include "../mtx.hpp"
#include "../grd.hpp"
#include "../vel.hpp"

#include "sdm_enums.hpp"

#include <memory>
using std::unique_ptr;

#include <map>
using std::map;

#include <boost/ptr_container/ptr_vector.hpp>
using boost::ptr_vector;

#include <boost/ptr_container/ptr_unordered_map.hpp>
using boost::ptr_unordered_map;

#include <vector>
using std::vector;

#include <string>
using std::string;

// TODO: option to set the number of threads to use!

/// @brief 
/// implementation of the Super-Droplet Method (@copydetails Shima_et_al_2009, QJRMS, 135)
/// with kappa-Koehler parameterisation of aerosol solubility (@copydetails Petters_and_Kreidenweis_2007, ACP, 7)
/// and ...

namespace sdm
{

  template <typename real_t>
  class sdm_proto
  {
    public: virtual void adjustments(
      mtx::arr<real_t> &rhod_th,
      mtx::arr<real_t> &rhod_rv,
      ptr_unordered_map<string, mtx::arr<real_t>> &aux, 
      const ptr_vector<mtx::arr<real_t>> C,
      const quantity<si::time, real_t> dt,
      bool record
    )  
    {
      assert(false);
    }
  };

  template <typename real_t, int thrust_device_system>
  class sdm : public sdm_proto<real_t>
  {
    private: struct detail;
    private: unique_ptr<detail> pimpl;

    // ctor 
    public: sdm(
      const grd<real_t> &grid, 
      const vel<real_t> &velocity,
      map<enum processes, bool> opts,
      enum xi_dfntns xi_dfntn,
      enum ode_algos adve_algo, enum ode_algos sedi_algo, enum ode_algos cond_algo, enum ode_algos chem_algo,
      const int adve_sstp, const int sedi_sstp, const int cond_sstp, const int chem_sstp, const int coal_sstp,
      real_t sd_conc_mean,
      real_t min_rd,
      real_t max_rd, 
      real_t mean_rd1, // dry aerosol initial distribution parameters
      real_t mean_rd2, // TODO: encapsulate in a map/enum
      real_t sdev_rd1,
      real_t sdev_rd2,
      real_t n1_tot,
      real_t n2_tot,
      real_t kappa,
      //initial chemical conditions (in air and droplets)
      map<enum chem_gas, quantity<phc::mixing_ratio, real_t>> opt_gas, 
      map<enum chem_aq, quantity<si::mass, real_t>> opt_aq 
    );

    public: void adjustments(
      mtx::arr<real_t> &rhod_th,
      mtx::arr<real_t> &rhod_rv,
      ptr_unordered_map<string, mtx::arr<real_t>> &aux, 
      const ptr_vector<mtx::arr<real_t>> C,
      const quantity<si::time, real_t> dt,
      bool record
    );
  };
};
