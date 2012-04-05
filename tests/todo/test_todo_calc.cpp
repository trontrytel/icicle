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

#include "../../src/phc.hpp"

typedef float real_t;

// simulation parameteters (the 8th WMO Cloud Modelling Workshop: Case 1    
// by W.W.Grabowski: http://rap.ucar.edu/~gthompsn/workshop2012/case1/case1.pdf  
const size_t 
  nx = 150,                   // 75 
  ny = 75;                    // 75 
const quantity<si::length,real_t> 
  dx = 20 * si::metres,       // 20 m
  dy = 20 * si::metres,       // 20 m
  z_inv = real_t(ny) * dy;
const quantity<si::temperature, real_t>
  th_0 = 289 * si::kelvins;   // 289 K
const quantity<phc::mixing_ratio, real_t>
  rt_0 = 7.5e-3;              // 7.5e-3 kg/kg
const quantity<si::pressure, real_t> 
  p_0 = 101500 * si::pascals; // 1015 hPa

// other parameters:
const int 
  bits = 32,
  fct = 1,  
  iord = 2;  
const size_t
  nout = 50,
  nt = 100;
const quantity<si::time, real_t> 
  dt = 10 * si::seconds; 

// pressure profile derived by integrating the hydrostatic eq.
// assuming constant theta, constant rv and R=R(rv) 
quantity<si::pressure,real_t> p_hydrostatic(
  quantity<si::length, real_t> z,
  quantity<si::temperature,real_t> th_0,
  quantity<phc::mixing_ratio, real_t> r,
  quantity<si::length, real_t> z_0,
  quantity<si::pressure, real_t> p_0
)
{
  return phc::p_1000<real_t>() * real_t(pow(
    pow(p_0 / phc::p_1000<real_t>(), phc::R_d_over_c_pd<real_t>()) 
    - 
    phc::R_d_over_c_pd<real_t>() * phc::g<real_t>() / th_0 / phc::R<real_t>(r) * (z - z_0), 
    phc::c_pd<real_t>() / phc::R_d<real_t>()
  ));
}

int main()
{
  {
    notice_macro("creating ini.nc ...")
    NcFile nf("ini.nc", NcFile::newFile);

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
      quantity<si::mass_density, real_t> rhod = (p - phc::p_v<real_t>(p, rt_0)) / (phc::exner<real_t>(p, rt_0) * phc::R_d<real_t>() * th_0);
      quantity<si::temperature, real_t> T = p / rhod / phc::R_d<real_t>();
      // theta^\star as a function of theta
      quantity<si::temperature, real_t> th = real_t(pow(
        real_t(th_0 / T) * pow(T / si::kelvins, pow(phc::R_over_c_p<real_t>(rt_0),-1)), 
        phc::R_over_c_p<real_t>(rt_0)
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
    << " --eqs todo"
    << " --grd.dx " << dx / si::metres
    << " --grd.nx " << nx
    << " --grd.dy " << dy / si::metres
    << " --grd.ny " << ny
    << " --adv mpdata"
      << " --adv.mpdata.fct " << fct
      << " --adv.mpdata.iord " << iord
    << " --vel rasinski"
      << " --vel.rasinski.Z_clb " << real_t(.5) * z_inv / si::metres
      << " --vel.rasinski.Z_top " << z_inv
      << " --vel.rasinski.A " << 1000
    << " --nt " << nt 
    << " --dt " << dt / si::seconds 
    << " --nout " << nout
    << " --out netcdf" 
    << " --out.netcdf.file out.nc"
    << " --slv openmp --nsd 2";
  system(cmd.str().c_str());
}
