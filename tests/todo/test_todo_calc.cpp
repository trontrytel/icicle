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
  nt = 10000;
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
    NcFile nf("ini.nc", NcFile::newFile);

    // dimensions
    NcDim ndx = nf.addDim("X", nx); 
    NcDim ndy = nf.addDim("Y", ny); 
    NcDim ndz = nf.addDim("Z", 1); // TODO: should not be needed
 
    // advected fields
    NcVar nvrhod_rv = nf.addVar("rhod_rv", ncFloat, vector<NcDim>({ndx, ndy}));
    NcVar nvrhod_rl = nf.addVar("rhod_rl", ncFloat, vector<NcDim>({ndx, ndy}));
    NcVar nvrhod_th = nf.addVar("rhod_th", ncFloat, vector<NcDim>({ndx, ndy}));

    // auxiliary fields
    NcVar nvrhod = nf.addVar("rhod", ncFloat, vector<NcDim>({ndx, ndy}));

    // initial values: (assuming no liquid water -> to be adjusted by the model)
    Array<real_t, 1> rhod(ny);
    for (int j = 0; j < ny; ++j) 
    {
      quantity<si::pressure, real_t> p =
        p_hydrostatic(real_t(j + .5) * si::metres, th_0, rt_0, real_t(0) * si::metres, p_0);
      rhod(j) = (
        (p - phc::p_v<real_t>(p, rt_0)) / (phc::exner<real_t>(p, rt_0) * phc::R_d<real_t>() * th_0) 
      ) / si::kilograms * si::metres * si::metres * si::metres;
    }
    Range j(0, ny-1);
    Array<real_t, 1> rhod_rl(j); 
    rhod_rl = 0; // to be adjusted by the model
    Array<real_t, 1> rhod_rv(rhod(j) * (rt_0 - rhod_rl(j)));
    Array<real_t, 1> rhod_th(rhod(j) * (th_0 / si::kelvins)); // TODO! theta^\star

cerr << "rhod: " << rhod << endl;

cerr << "rhod_rv: " << rhod_rv << endl;

cerr << "rhod_th: " << rhod_th << endl;

    // writing the profiles to the netCDF 
    for (size_t i = 0; i < nx; ++i) // TODO: this logic should be handled within icicle! - it's quite usual to provide profiles
    {
      nvrhod.putVar(start({i,0}), count({1,ny}), rhod.data());
      nvrhod_rv.putVar(start({i,0}), count({1,ny}), rhod_rv.data());
      nvrhod_rl.putVar(start({i,0}), count({1,ny}), rhod_rl.data());
      nvrhod_th.putVar(start({i,0}), count({1,ny}), rhod_th.data());
    }
  }
  
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
