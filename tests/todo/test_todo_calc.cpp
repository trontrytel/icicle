/* @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
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

#include <iostream>
using std::ostringstream;
using std::cerr;
using std::endl;
#define notice_macro(msg) { cerr << msg << endl; }

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

#include "../../src/phc.hpp"

typedef float real_t;

// simulation parameteters (the 8th WMO Cloud Modelling Workshop: Case 1    
// by W.W.Grabowski: http://rap.ucar.edu/~gthompsn/workshop2012/case1/case1.pdf  
const size_t 
  nx = 75,                   // 75 
  ny = 75;                    // 75 
const quantity<si::length,real_t> 
  dx = 20 * si::metres,       // 20 m
  dy = 20 * si::metres;       // 20 m
const quantity<si::temperature, real_t>
  th_0 = 289 * si::kelvins;   // 289 K
const quantity<phc::mixing_ratio, real_t>
  rt_0 = 7.5e-3;              // 7.5e-3 kg/kg
const quantity<si::pressure, real_t> 
  p_0 = 101500 * si::pascals; // 1015 hPa

// other parameters deduced from the Fortran code published at:
// http://www.rap.ucar.edu/~gthompsn/workshop2012/case1/kinematic_wrain.vocals.v3.for
const int 
  bits = 64,
  fct = 0,  
  iord = 2;  
const quantity<si::time, real_t> 
  t_max = 3600 * si::seconds, // 4 * 3600
  dt_out = real_t(50) * si::seconds; // 300
const quantity<si::velocity, real_t>
  w_max = real_t(.6) * si::metres / si::second; // TODO: check it!
const quantity<si::mass_density, real_t>
  rho_0 = 1 * si::kilograms / si::cubic_metres;
const quantity<divide_typeof_helper<si::momentum, si::area>::type, real_t> 
  ampl = rho_0 * w_max * (real_t(nx) * dx) / real_t(4*atan(1));

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

int main()
{
  {
    notice_macro("creating rho.nc ...")
    NcFile nf("rho.nc", NcFile::newFile, NcFile::classic);

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
  }
  {
    notice_macro("creating ini.nc ...")
    NcFile nf("ini.nc", NcFile::newFile, NcFile::classic);

    // dimensions
    NcDim ndx = nf.addDim("X", 1); 
    NcDim ndy = nf.addDim("Y", ny); 
    NcDim ndz = nf.addDim("Z", 1); 

    // dimension-annotating variable
    NcVar nvy = nf.addVar("Y", ncFloat, vector<NcDim>({ndx,ndy,ndz}));
 
    // advected fields
    NcVar nvrhod_rv = nf.addVar("rhod_rv", ncFloat, vector<NcDim>({ndx,ndy,ndz}));
    NcVar nvrhod_rl = nf.addVar("rhod_rl", ncFloat, vector<NcDim>({ndx,ndy,ndz}));
    NcVar nvrhod_th = nf.addVar("rhod_th", ncFloat, vector<NcDim>({ndx,ndy,ndz}));

    // auxiliary fields
    NcVar nvrhod = nf.addVar("rhod", ncFloat, vector<NcDim>({ndx,ndy,ndz}));

    // initial values: (assuming no liquid water -> to be adjusted by the model)
    Array<real_t, 1> arr_rhod(ny), arr_z(ny), arr_rhod_rl(ny), arr_rhod_rv(ny), arr_rhod_th(ny);
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
      arr_rhod_rv(j) = arr_rhod(j) * rt_0;
      arr_rhod_th(j) = arr_rhod(j) * th / si::kelvins;
    }

    // writing the profiles to the netCDF 
    nvy.putVar(arr_z.data());
    nvrhod.putVar(arr_rhod.data());
    nvrhod_rv.putVar(arr_rhod_rv.data());
    nvrhod_rl.putVar(arr_rhod_rl.data());
    nvrhod_th.putVar(arr_rhod_th.data());
  }
  
  notice_macro("calling the solver ...")
  ostringstream cmd;
  cmd
    << "../../icicle"
    << " --bits " << bits
    << " --ini netcdf"
    << " --ini.netcdf.file ini.nc"
    << " --eqs todo_bulk"
    << " --grd.dx " << dx / si::metres
    << " --grd.nx " << nx
    << " --grd.dy " << dy / si::metres
    << " --grd.ny " << ny
    << " --adv mpdata"
      << " --adv.mpdata.fct " << fct
      << " --adv.mpdata.iord " << iord
    << " --vel rasinski"
      << " --vel.rasinski.file " << "rho.nc"
      << " --vel.rasinski.A " << ampl
    << " --t_max " << real_t(t_max / si::seconds)
    << " --dt " << "auto" 
    << " --dt_out " << real_t(dt_out / si::seconds)
    << " --out netcdf" 
    << " --out.netcdf.file out.nc"
    //<< " --slv serial"
    << " --slv openmp --nsd 3"
    ;
  system(cmd.str().c_str());
}
