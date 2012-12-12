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

#include <cgicc/Cgicc.h>
cgicc::Cgicc cgi;

template <typename T>
T http_or_default(const string &name, const T &def)
{
  static char *qs = getenv("QUERY_STRING");
  if (qs == NULL) return def;
  string tmp;
  return boost::lexical_cast<T>(cgi(name));
}

typedef float real_t;

// simulation parameteters (the 8th WMO Cloud Modelling Workshop: Case 1    
// by W.W.Grabowski: http://rap.ucar.edu/~gthompsn/workshop2012/case1/case1.pdf  
const size_t 
  nx = http_or_default("nx",size_t(76)),                   // 75 
  ny = http_or_default("ny",size_t(76));                    // 75 
const quantity<si::length,real_t> 
  dx = http_or_default("dx",real_t(20)) * si::metres,       // 20 m
  dy = http_or_default("dy",real_t(20)) * si::metres;       // 20 m
const quantity<si::temperature, real_t>
  th_0 = http_or_default("th_0", real_t(289)) * si::kelvins;   // 289 K
const quantity<phc::mixing_ratio, real_t>
  rt_0 = http_or_default("rt_0", real_t(7.5e-3));              // 7.5e-3 kg/kg

// other parameters deduced from the Fortran code published at:
// http://www.rap.ucar.edu/~gthompsn/workshop2012/case1/kinematic_wrain.vocals.v3.for
const int 
  bits = http_or_default("bits", int(64)),
  fct = http_or_default("fct", int(0)),  
  toa = http_or_default("toa", int(0)),
  iord = http_or_default("iord", int(2)),
  nsd = http_or_default("nsd", int(4));  
const quantity<si::time, real_t> 
  t_max = 50 * si::seconds, // 4 * 3600
  dt_out = real_t(5) * si::seconds; // 300
const quantity<si::velocity, real_t>
  w_max = http_or_default("w_max", real_t(.6)) * si::metres_per_second; // .6 TODO: check it!

// options for microphysics
std::string micro_opt = http_or_default("micro", string("sdm")); 

// blk parameters
bool 
  blk_cevp = http_or_default("blk_cevp", true),
  blk_conv = http_or_default("blk_conv", true),
  blk_clct = http_or_default("blk_clct", true),
  blk_sedi = http_or_default("blk_sedi", true),
  blk_revp = http_or_default("blk_revp", true);

// mm parameters
bool
  mm_act   = true,
  mm_cond  = true,
  mm_autoc = true,
  mm_acc   = true,
//  mm_self  = false, no selfaccretion in KK2000
  mm_turb  = false,
  mm_sedi  = false;
real_t
  mean_rd = .04e-6,  //TODO - bimodal aerosol size spectrum
  sdev_rd = 1.4,
  n_tot = 60e6;

int
  sdm_adve_sstp = 1,
  sdm_sedi_sstp = 1,
  sdm_chem_sstp = 100,
  sdm_cond_sstp = 1,
  sdm_coal_sstp = 1;
bool 
  sdm_adve = http_or_default("sdm_adve", true),
  sdm_cond = http_or_default("sdm_cond", true),
  sdm_coal = http_or_default("sdm_coal", true),
  sdm_sedi = http_or_default("sdm_sedi", true),
  sdm_chem = http_or_default("sdm_chem", true);
real_t 
  sd_conc_mean = http_or_default("sd_conc_mean", 256);

using wkc::icmw8_case1::micro_t;
using wkc::icmw8_case1::bulk;
using wkc::icmw8_case1::sdm;
using wkc::icmw8_case1::mm;

int main(int argc, char **argv)
{
  enum micro_t micro;
  if (micro_opt == "bulk") micro = bulk;
  else if (micro_opt == "sdm") micro = sdm;
  else if (micro_opt == "mm") micro = mm;
  else error_macro("unknown microphysics")

  string dir = argc > 1 ? argv[1] : "tmp";
  notice_macro("creating " << dir << "...");
  boost::filesystem::create_directory(dir);

  notice_macro("creating input files ...")
  string opts = wkc::icmw8_case1::create_files(
    micro,
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

    << " --t_max " << real_t(t_max / si::seconds)
    << " --dt " << "auto" 
    << " --dt_out " << real_t(dt_out / si::seconds)

// TODO TEMP TODO TEMP for mm !!!
//    << " --dt " << real_t(1.)
//    << " --nt " << real_t(3000)  //4*3600
//    << " --nout " << real_t(100)   //600

    << " --out netcdf" 
    << " --out.netcdf.file " << dir << "/out.nc";

    // TEMP TODO!
    if (micro == bulk) 
      cmd << " --slv openmp --nsd " << nsd;
    else 
      cmd << " --slv serial";

    if (micro == bulk) cmd 
      << " --eqs.todo_bulk.cevp " << blk_cevp
      << " --eqs.todo_bulk.conv " << blk_conv
      << " --eqs.todo_bulk.clct " << blk_clct
      << " --eqs.todo_bulk.sedi " << blk_sedi
      << " --eqs.todo_bulk.revp " << blk_revp
    ;
    else if (micro == mm) cmd 
      << " --eqs.todo_mm.act "   << mm_act
      << " --eqs.todo_mm.cond "  << mm_cond
      << " --eqs.todo_mm.acc "   << mm_acc
      << " --eqs.todo_mm.autoc " << mm_autoc
      << " --eqs.todo_mm.turb "  << mm_turb
      << " --eqs.todo_mm.sedi "  << mm_sedi
    ;
    else if (micro == sdm) cmd 
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
    ;
    else assert(false);
    
  if (EXIT_SUCCESS != system(cmd.str().c_str()))
    error_macro("model run failed: " << cmd.str())
}
