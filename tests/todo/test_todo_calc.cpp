/* @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date March 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#include <netcdf>
using netCDF::NcFile;
using netCDF::NcDim;
using netCDF::NcVar;
using netCDF::ncFloat;

#include <vector>
using std::vector;

#include <iostream>
using std::ostringstream;

#include <blitz/array.h>
using blitz::Array;
using blitz::Range;
using blitz::toEnd;

#include <boost/units/systems/si.hpp>
#include <boost/units/io.hpp>
namespace si = boost::units::si;
using boost::units::quantity;

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
  th_l_0 = 289 * si::kelvins; // 289 K
const quantity<si::dimensionless, real_t>
  q_t_0 = 7.5e-3;             // 7.5e-3 kg/kg

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

    // initial condition 
    Array<real_t, 2> rhod_l(nx, ny);
    rhod_l = 0;
    rhod_l(Range::all(), Range(ny/2, toEnd)) = .001;
    nvrhod_rl.putVar(rhod_l.data());
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
