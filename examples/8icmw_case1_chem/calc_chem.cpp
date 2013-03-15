/* @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @author Ania Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date March 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#include <string>
using std::string;

#include <iostream>
using std::ostringstream;
using std::cerr;
using std::endl;

#include <boost/units/systems/si.hpp>
#include <boost/units/io.hpp>
namespace si = boost::units::si;
using boost::units::quantity;
using boost::units::divide_typeof_helper;
using boost::units::detail::get_value;

#include <boost/filesystem.hpp>

#include <boost/lexical_cast.hpp>

#include "../../src/cmn/cmn_error.hpp"
#include "../../src/cmn/cmn_units.hpp"
#include "../../src/phc/phc.hpp"
#include "../../src/wkc/wkc_icmw8_case1.hpp"

typedef float real_t;

// simulation parameteters (the 8th WMO Cloud Modelling Workshop: Case 1    
// by W.W.Grabowski: http://rap.ucar.edu/~gthompsn/workshop2012/case1/case1.pdf  
const size_t 
  nx = size_t(76),                   // 75 
  ny = size_t(76);                    // 75 
const quantity<si::length,real_t> 
  dx = real_t(20) * si::metres,       // 20 m
  dy = real_t(20) * si::metres;       // 20 m
const quantity<si::temperature, real_t>
  th_0 = real_t(289) * si::kelvins;   // 289 K
const quantity<phc::mixing_ratio, real_t>
  rt_0 =real_t(7.5e-3);              // 7.5e-3 kg/kg

// other parameters deduced from the Fortran code published at:
// http://www.rap.ucar.edu/~gthompsn/workshop2012/case1/kinematic_wrain.vocals.v3.for
const int 
  bits = int(64),
  fct = int(0),  
  toa = int(0),
  iord = int(2);
const quantity<si::velocity, real_t>
  w_max = real_t(.6) * si::metres_per_second; // .6 TODO: check it!

// options for microphysics
std::string micro_opt = string("sdm"); 

int
  sdm_adve_sstp = 1,
  sdm_sedi_sstp = 1,
  sdm_chem_sstp = 5,
  sdm_cond_sstp = 5,
  sdm_coal_sstp = 1;
bool 
  sdm_adve = true,
  sdm_cond = true,
  sdm_coal = true,
  sdm_sedi = true,
  sdm_chem = true;
real_t 
  sd_conc_mean = 128;

using wkc::icmw8_case1::micro_t;
using wkc::icmw8_case1::sdm;

int main(int argc, char **argv)
{
  string dir = argc > 1 ? argv[1] : "tmp";
  notice_macro("creating " << dir << "...");
  boost::filesystem::create_directory(dir);

  notice_macro("creating input files ...")
  string opts = wkc::icmw8_case1::create_files(
    sdm,
    dir + "/rho.nc", 
    dir + "/ini.nc",
    nx, ny, dx, dy, th_0, rt_0, w_max
  );

  notice_macro("calling the solver ...")
  ostringstream cmd;
  cmd
    << "../../icicle " << opts 
    << " --adv mpdata"
      << " --adv.mpdata.fct " << fct
      << " --adv.mpdata.iord " << iord
      << " --adv.mpdata.third_order " << toa

    << " --dt " << real_t(1)
    << " --nt " << real_t(3600) // 1800
    << " --nout " << real_t(30) // 60

    << " --out netcdf" 
    << " --out.netcdf.file " << dir << "/out.nc"

    << " --slv serial"

    << " --eqs.todo_sdm.adve.sstp " << sdm_adve_sstp
    << " --eqs.todo_sdm.sedi.sstp " << sdm_sedi_sstp
    << " --eqs.todo_sdm.cond.sstp " << sdm_cond_sstp
    << " --eqs.todo_sdm.chem.sstp " << sdm_chem_sstp
    << " --eqs.todo_sdm.coal.sstp " << sdm_coal_sstp

    << " --eqs.todo_sdm.adve " << sdm_adve
    << " --eqs.todo_sdm.cond " << sdm_cond
    << " --eqs.todo_sdm.coal " << sdm_coal
    << " --eqs.todo_sdm.sedi " << sdm_sedi
    << " --eqs.todo_sdm.chem " << sdm_chem
    << " --eqs.todo_sdm.sd_conc_mean " << sd_conc_mean

    // TEMP...
    << " --eqs.todo_sdm.out_m0 \"";
      for (int i = 0; i < 25; ++i) cmd << i << "e-6:" << i + 1 << "e-6;";
      for (int i = 0; i < 7; ++i) cmd << 25 * int(pow(2,i)) << "e-6:" << 25 * int(pow(2,i+1)) << "e-6;";
      cmd << "\"";
    
  if (EXIT_SUCCESS != system(cmd.str().c_str()))
    error_macro("model run failed: " << cmd.str())
}
