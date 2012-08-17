/* @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @author Ania Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date March 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#include <vector>
using std::vector;

#include <netcdf>
using netCDF::NcFile;
using netCDF::NcDim;
using netCDF::NcVar;
using netCDF::ncFloat;
typedef vector<size_t> start;
typedef vector<size_t> count;

#include <string>
using std::string;

#include <iostream>
using std::ostringstream;
using std::cerr;
using std::endl;

#include <blitz/array.h>
using blitz::Array;
using blitz::Range;
using blitz::toEnd;

#include <boost/units/systems/si.hpp>
#include <boost/units/io.hpp>
namespace si = boost::units::si;
using boost::units::quantity;
using boost::units::divide_typeof_helper;
using boost::units::detail::get_value;

// mkdir()
#include <sys/stat.h>
#include <sys/types.h>

#include <boost/lexical_cast.hpp>

#include "../../src/cmn.hpp"
#include "../../src/phc_theta.hpp"
#include "../../src/phc_terminal_vel.hpp"

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
  nx = http_or_default("nx",size_t(75)),                   // 75 
  ny = http_or_default("ny",size_t(75));                    // 75 
const quantity<si::length,real_t> 
  dx = http_or_default("dx",real_t(20.)) * si::metres,       // 20 m
  dy = http_or_default("dy",real_t(20.)) * si::metres;       // 20 m
const quantity<si::temperature, real_t>
  th_0 = http_or_default("th_0", real_t(289)) * si::kelvins;   // 289 K
const quantity<phc::mixing_ratio, real_t>
  rt_0 = http_or_default("rt_0", real_t(7.5e-3));              // 7.5e-3 kg/kg
const quantity<si::pressure, real_t> 
  p_0 = 101500 * si::pascals; // 1015 hPa
const quantity<si::dimensionless, real_t>
  kappa = 0.61; // ammonium sulphate
  //kappa = 1.28; // sodium chloride
const quantity<si::dimensionless, real_t>// equivalent of kappa for 2 moment param. (assuming no insoluable aerosol core)
  //TODO? calculate from kappa?
  chem_b = 0.505; // ammonium sluphate
  //chem_b = 1.33; // sodium chloride

// other parameters deduced from the Fortran code published at:
// http://www.rap.ucar.edu/~gthompsn/workshop2012/case1/kinematic_wrain.vocals.v3.for
const int 
  bits = http_or_default("bits", int(32)),
  fct = http_or_default("fct", int(0)),  
  toa = http_or_default("toa", int(0)),
  iord = http_or_default("iord", int(1)),
  nsd = http_or_default("nsd", int(1));  
const quantity<si::time, real_t> 
  t_max = 400 * si::seconds, // 4 * 3600
  dt_out = real_t(3) * si::seconds; // 300
const quantity<si::velocity, real_t>
  w_max = http_or_default("w_max", real_t(.6)) * si::metres / si::second; // .6 TODO: check it!
const quantity<si::mass_density, real_t>
  rho_0 = 1 * si::kilograms / si::cubic_metres;
const quantity<divide_typeof_helper<si::momentum, si::area>::type, real_t> 
  ampl = rho_0 * w_max * (real_t(nx) * dx) / real_t(4*atan(1));

// options for microphysics
std::string micro = http_or_default("micro", string("mm")); // sdm | bulk | mm

// blk parameters
bool 
  blk_cevp = http_or_default("blk_cevp", true),
  blk_conv = http_or_default("blk_conv", true),
  blk_clct = http_or_default("blk_clct", true),
  blk_sedi = http_or_default("blk_sedi", true),
  blk_revp = http_or_default("blk_revp", true);

// mm parameters
bool
  mm_act   = "true",
  mm_cond  = "false",
  mm_acc   = "false",
  mm_autoc = "false",
  mm_self  = "false",
  mm_turb  = "false",
  mm_sedi  = "false";
real_t
  mean_rd = .04e-6,  //TODO - bimodal aerosol size spectrum
  sdev_rd = 1.4,
  n_tot = 60e6;

// sdm parameters
std::string
  sdm_xi = "p2", 
  sdm_ode_algo_adve = "euler",
  sdm_ode_algo_sedi = "euler",
  sdm_ode_algo_cond = "euler",
  sdm_ode_algo_chem = "euler";
bool 
  sdm_adve = http_or_default("sdm_adve", true),
  sdm_cond = http_or_default("sdm_cond", true),
  sdm_coal = http_or_default("sdm_coal", true),
  sdm_sedi = http_or_default("sdm_sedi", true),
  sdm_chem = http_or_default("sdm_chem", true);
real_t 
  sd_conc_mean = http_or_default("sd_conc_mean", 32),
  mean_rd1 = .04e-6,
  mean_rd2 = .15e-6,
  sdev_rd1 = 1.4,
  sdev_rd2 = 1.6,
  n1_tot = 60e6,
  n2_tot = 40e6,
  min_rd = 1e-9,
  max_rd = 1e-6;

// pressure profile derived by integrating the hydrostatic eq.
// assuming constant theta, constant rv and R=R(rv) 
quantity<si::pressure,real_t> p_hydrostatic(
  quantity<si::length, real_t> z,
  quantity<si::temperature,real_t> th_0,
  quantity<phc::mixing_ratio, real_t> r_0,
  quantity<si::length, real_t> z_0,
  quantity<si::pressure, real_t> p_0
)
{
  return phc::p_1000<real_t>() * real_t(pow(
    pow(p_0 / phc::p_1000<real_t>(), phc::R_d_over_c_pd<real_t>()) 
    - 
    phc::R_d_over_c_pd<real_t>() * phc::g<real_t>() / th_0 / phc::R<real_t>(r_0) * (z - z_0), 
    phc::c_pd<real_t>() / phc::R_d<real_t>()
  ));
}

// TODO: document it and name it correctly
quantity<si::mass_density, real_t> rhod_todo(
  quantity<si::pressure, real_t> p,
  quantity<si::temperature, real_t> th_0,
  quantity<phc::mixing_ratio, real_t> r_0
)
{
  return (p - phc::p_v<real_t>(p, rt_0)) / (phc::exner<real_t>(p, rt_0) * phc::R_d<real_t>() * th_0);
}

int main(int argc, char **argv)
{
  string dir = argc > 1 ? argv[1] : "tmp";
  notice_macro("creating " << dir << "...");
  mkdir(dir.c_str(), 0777);

  {
    notice_macro("creating rho.nc ...")
    NcFile nf(dir + "/rho.nc", NcFile::newFile, NcFile::classic);

    // dimensions and variables
    NcDim ndy = nf.addDim("Y", 2 * ny+1); 
    NcVar nvy = nf.addVar("Y", ncFloat, vector<NcDim>({ndy}));
    NcVar nvrho = nf.addVar("rho", ncFloat, vector<NcDim>({ndy}));
  
    // dry air density at all needed levels
    for (size_t j = 0; j < 2 * ny + 1; ++j) 
    {
      // calculating
      quantity<si::length, real_t> z = real_t(.5 * j) * dy ;
      quantity<si::pressure, real_t> p = p_hydrostatic(z, th_0, rt_0, real_t(0) * si::metres, p_0);
      quantity<si::mass_density, real_t> rhod = rhod_todo(p, th_0, rt_0);

      // writing to the netCDF
      nvy.putVar(start({j}), z / si::metres);
      nvrho.putVar(start({j}), rhod * si::cubic_metres / si::kilogram);
    }
  } // nf gets closed here

  {
    notice_macro("creating ini.nc ...")
    NcFile nf(dir + "/ini.nc", NcFile::newFile, NcFile::classic);

    // dimensions
    NcDim ndx = nf.addDim("X", 1); 
    NcDim ndy = nf.addDim("Y", ny); 
    NcDim ndz = nf.addDim("Z", 1); 

    // dimension-annotating variable
    NcVar nvy = nf.addVar("Y", ncFloat, vector<NcDim>({ndx,ndy,ndz}));
 
    // advected fields
    NcVar nvrhod_rv = nf.addVar("rhod_rv", ncFloat, vector<NcDim>({ndx,ndy,ndz}));
    NcVar nvrhod_rl = nf.addVar("rhod_rl", ncFloat, vector<NcDim>({ndx,ndy,ndz}));
    NcVar nvrhod_rr = nf.addVar("rhod_rr", ncFloat, vector<NcDim>({ndx,ndy,ndz}));
    NcVar nvrhod_th = nf.addVar("rhod_th", ncFloat, vector<NcDim>({ndx,ndy,ndz}));
    
    NcVar nvrhod_nl, nvrhod_nr;
    if (micro == "mm")
    {
      nvrhod_nl = nf.addVar("rhod_nl", ncFloat, vector<NcDim>({ndx,ndy,ndz}));
      nvrhod_nr = nf.addVar("rhod_nr", ncFloat, vector<NcDim>({ndx,ndy,ndz}));
    }
    // auxiliary fields
    NcVar nvrhod = nf.addVar("rhod", ncFloat, vector<NcDim>({ndx,ndy,ndz}));

    // initial values: (assuming no liquid and rain water -> to be adjusted by the model)
    Array<real_t, 1> arr_rhod(ny), arr_z(ny), arr_rhod_rl(ny), 
                     arr_rhod_rv(ny), arr_rhod_rr(ny), arr_rhod_th(ny),
                     arr_rhod_nl(ny), arr_rhod_nr(ny);
    for (int j = 0; j < ny; ++j) 
    {
      quantity<si::length, real_t> z = real_t(j + .5) * dy;
      quantity<si::pressure, real_t> p = p_hydrostatic(z, th_0, rt_0, real_t(0) * si::metres, p_0);
      quantity<si::mass_density, real_t> rhod = rhod_todo(p, th_0, rt_0); 
      quantity<si::temperature, real_t> T = p / rhod / phc::R_d<real_t>();
      // theta^\star as a function of theta
      quantity<si::temperature, real_t> th = real_t(pow(
        real_t(th_0 / T) * pow(T / si::kelvins, pow(phc::R_over_c_p<real_t>(rt_0) / phc::R_d_over_c_pd<real_t>(),-1)), 
        phc::R_over_c_p<real_t>(rt_0) / phc::R_d_over_c_pd<real_t>()
      )) * si::kelvins;

      arr_z(j) = z / si::metres;
      arr_rhod(j) = rhod / si::kilograms * si::cubic_metres;
      arr_rhod_rl(j) = 0; // to be adjusted by the model
      arr_rhod_rr(j) = 0; // to be adjusted by the model
      arr_rhod_rv(j) = arr_rhod(j) * rt_0;
      arr_rhod_th(j) = arr_rhod(j) * th / si::kelvins;
      if (micro == "mm"){
        arr_rhod_nl(j)=0; // to be adjusted by the model
        arr_rhod_nr(j)=0; // to be adjusted by the model
      }    
    }

    // writing the profiles to the netCDF 
    nvy.putVar(arr_z.data());
    nvrhod.putVar(arr_rhod.data());
    nvrhod_rv.putVar(arr_rhod_rv.data());
    nvrhod_rl.putVar(arr_rhod_rl.data());
    nvrhod_rr.putVar(arr_rhod_rr.data());
    nvrhod_th.putVar(arr_rhod_th.data());
    if (micro == "mm"){
      nvrhod_nl.putVar(arr_rhod_nl.data());
      nvrhod_nr.putVar(arr_rhod_nr.data());
    }
  }
  
  notice_macro("calling the solver ...")
  ostringstream cmd;
  cmd
    << "../../icicle"
    << " --bits " << bits
    << " --ini netcdf"
    << " --ini.netcdf.file " << dir << "/ini.nc"
    << " --grd.dx " << dx / si::metres
    << " --grd.nx " << nx
    << " --grd.dy " << dy / si::metres
    << " --grd.ny " << ny
    << " --adv mpdata"
      << " --adv.mpdata.fct " << fct
      << " --adv.mpdata.iord " << iord
      << " --adv.mpdata.third_order " << toa
    << " --vel rasinski"
      << " --vel.rasinski.file " << dir << "/rho.nc"
      << " --vel.rasinski.A " << ampl
    << " --t_max " << real_t(t_max / si::seconds)
    << " --dt " << "auto" 
    << " --dt_out " << real_t(dt_out / si::seconds)
    << " --out netcdf" 
    << " --out.netcdf.file " << dir << "/out.nc";
    if (micro == "bulk") 
      cmd << " --slv openmp --nsd " << nsd;
    else 
      cmd << " --slv serial";
    if (micro == "bulk") cmd << " --eqs todo_bulk"
      << " --eqs.todo_bulk.cevp " << blk_cevp
      << " --eqs.todo_bulk.conv " << blk_conv
      << " --eqs.todo_bulk.clct " << blk_clct
      << " --eqs.todo_bulk.sedi " << blk_sedi
      << " --eqs.todo_bulk.revp " << blk_revp
    ;
    else if (micro == "mm") cmd << " --eqs todo_mm"
      << " --eqs.todo_mm.act "   << mm_act
      << " --eqs.todo_mm.cond "  << mm_cond
      << " --eqs.todo_mm.acc "   << mm_acc
      << " --eqs.todo_mm.autoc " << mm_autoc
      << " --eqs.todo_mm.self "  << mm_self
      << " --eqs.todo_mm.turb "  << mm_turb
      << " --eqs.todo_mm.sedi "  << mm_sedi
      << " --eqs.todo_mm.mean_rd " << mean_rd
      << " --eqs.todo_mm.sdev_rd " << sdev_rd
      << " --eqs.todo_mm.n_tot "   << n_tot
      << " --eqs.todo_mm.chem_b "    << chem_b
    ;
    else if (micro == "sdm") cmd << " --eqs todo_sdm"
      << " --eqs.todo_sdm.xi " << sdm_xi
      << " --eqs.todo_sdm.adve.algo " << sdm_ode_algo_adve
      << " --eqs.todo_sdm.sedi.algo " << sdm_ode_algo_sedi
      << " --eqs.todo_sdm.cond.algo " << sdm_ode_algo_cond
      << " --eqs.todo_sdm.chem.algo " << sdm_ode_algo_chem
      << " --eqs.todo_sdm.adve " << sdm_adve
      << " --eqs.todo_sdm.cond " << sdm_cond
      << " --eqs.todo_sdm.coal " << sdm_coal
      << " --eqs.todo_sdm.sedi " << sdm_sedi
      << " --eqs.todo_sdm.chem " << sdm_chem
      << " --eqs.todo_sdm.sd_conc_mean " << sd_conc_mean
      << " --eqs.todo_sdm.min_rd " << min_rd
      << " --eqs.todo_sdm.max_rd " << max_rd
      << " --eqs.todo_sdm.mean_rd1 " << mean_rd1
      << " --eqs.todo_sdm.mean_rd2 " << mean_rd2
      << " --eqs.todo_sdm.sdev_rd1 " << sdev_rd1
      << " --eqs.todo_sdm.sdev_rd2 " << sdev_rd2
      << " --eqs.todo_sdm.n1_tot " << n1_tot
      << " --eqs.todo_sdm.n2_tot " << n2_tot
      << " --eqs.todo_sdm.kappa " << kappa
    ;
    else assert(false);
    
  if (EXIT_SUCCESS != system(cmd.str().c_str()))
    cerr << "model run failed: " << cmd.str() << endl;
}
