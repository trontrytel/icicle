/* @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @author Ania Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date March 2013
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @section DESCRIPTION
 *    helper routines for the chem simulations (based on case 1 and Kreidenweis 2003)
 */

#pragma once

#include <string>
using std::string;

#include <iostream>
using std::ostringstream;
using std::cerr;
using std::endl;

#include <vector>
using std::vector;

#include <list>
using std::list;

#include <map>
using std::map;

// TODO: are all these "usings" neccessary
#include <netcdf>
using netCDF::NcFile;
using netCDF::NcDim;
using netCDF::NcVar;
using netCDF::ncFloat;
typedef vector<size_t> start;
typedef vector<size_t> count;

#include "../../src/phc/phc_theta.hpp"
//#include "../../src/phc/phc_terminal_vel.hpp"
#include "../../src/phc/phc_hydrostatic.hpp"

#include <blitz/array.h>
using blitz::Array;
using blitz::Range;
using blitz::toEnd;

#include "../cmn/cmn_error.hpp"

namespace wkc 
{
  namespace tmp_chem
  {
    // returns a string with options to be passed to icicle
    template <typename real_t>
    string create_files(
      const string rho_filename = "rho.nc",
      const string ini_filename = "ini.nc",
      const int nx = 75,
      const int ny = 75,
      const quantity<si::length, real_t> dx = real_t(20) * si::metres,
      const quantity<si::length, real_t> dy = real_t(20) * si::metres,
      const quantity<si::temperature, real_t> th_0 = real_t(289) * si::kelvins,
      const quantity<si::dimensionless, real_t> rt_0 = real_t(7.5e-3),
      const quantity<si::velocity, real_t> w_max = real_t(.6) * si::metres_per_second // TODO: check it
    ) 
    {
      ostringstream opts;

      opts 
        << " --grd.dx " << dx / si::metres
        << " --grd.nx " << nx
        << " --grd.dy " << dy / si::metres
        << " --grd.ny " << ny
        << " --eqs todo_sdm"
        // initial dry aerosol spectra as defined in case 1 of 8icmw
	<< " --eqs.todo_sdm.mean_rd1 " << .04e-6 / 2
	<< " --eqs.todo_sdm.mean_rd2 " << .15e-6 / 2    //not used (one mode only)
	<< " --eqs.todo_sdm.sdev_rd1 " << 1.4
	<< " --eqs.todo_sdm.sdev_rd2 " << 1.6           //not used (one mode only)
	<< " --eqs.todo_sdm.n1_tot " << 60e6
	<< " --eqs.todo_sdm.n2_tot " << 0               //not used (one mode only)
	<< " --eqs.todo_sdm.kappa " << 0.9; 
       //kappa = .61 // ammonium sulphate
       //kappa = 1.28; // sodium chloride NaCl
       //kappa = .9; // H2SO4

      // scope block introduced in order to close the netCDF file
      {
        notice_macro("creating " << rho_filename << " ...")
        NcFile nf(rho_filename, NcFile::newFile, NcFile::classic);

        // dimensions and variables
        NcDim ndy = nf.addDim("Y", 2 * ny+1);
        NcVar nvy = nf.addVar("Y", ncFloat, vector<NcDim>({ndy}));
        NcVar nvrho = nf.addVar("rho", ncFloat, vector<NcDim>({ndy}));

        // dry air density at all needed levels
        for (size_t j = 0; j < 2 * ny + 1; ++j)
        {
          // calculating
          // <TODO> repeated above!
          const quantity<si::pressure, real_t> p_0 = real_t(101500) * si::pascals;
          quantity<si::length, real_t> z = real_t(.5 * j) * dy ;
          quantity<si::pressure, real_t> p = phc::hydrostatic::p(z, th_0, rt_0, real_t(0) * si::metres, p_0);
          quantity<si::mass_density, real_t> rhod = phc::rhod(p, th_0, rt_0);
          // </TODO>

          // writing to the netCDF
          nvy.putVar(start({j}), z / si::metres);
          nvrho.putVar(start({j}), rhod * si::cubic_metres / si::kilogram);
        }
        notice_macro("done.")
      }

      {
        // TODO: check it!
        const quantity<si::mass_density, real_t>
          rho_0 = 1 * si::kilograms / si::cubic_metres;
        const quantity<divide_typeof_helper<si::momentum, si::area>::type, real_t>
          ampl = rho_0 * w_max * (real_t(nx) * dx) / phc::pi<real_t>();
 
        opts << " --vel rasinski"
          << " --vel.rasinski.file " << rho_filename
          << " --vel.rasinski.A " << ampl;
      }

      // scope block introduced in order to close the netCDF file
      {
        opts
          << " --ini netcdf"
          << " --ini.netcdf.file " << ini_filename;

	notice_macro("creating " << ini_filename << " ...")
	NcFile nf(ini_filename, NcFile::newFile, NcFile::classic);

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

	// auxiliary fields
	NcVar nvrhod = nf.addVar("rhod", ncFloat, vector<NcDim>({ndx,ndy,ndz}));

	// initial values: (assuming no liquid and rain water -> to be adjusted by the model)
	Array<real_t, 1> arr_rhod(ny), arr_z(ny), arr_rhod_rl(ny),
			 arr_rhod_rv(ny), arr_rhod_th(ny),
			 arr_zero(ny);
	for (int j = 0; j < ny; ++j)
	{
          // <TODO> repeated below!
          const quantity<si::pressure, real_t> p_0 = real_t(101500) * si::pascals;
	  quantity<si::length, real_t> z = real_t(j + .5) * dy;
	  quantity<si::pressure, real_t> p = phc::hydrostatic::p(z, th_0, rt_0, real_t(0) * si::metres, p_0);
	  quantity<si::mass_density, real_t> rhod = phc::rhod(p, th_0, rt_0);
          // </TODO>

	  quantity<si::temperature, real_t> T = p / rhod / phc::R_d<real_t>();
	  // theta^\star as a function of theta
	  quantity<si::temperature, real_t> th = real_t(pow(
	    real_t(th_0 / T) * pow(T / si::kelvins, pow(phc::R_over_c_p<real_t>(rt_0) / phc::R_d_over_c_pd<real_t>(),-1)),
	    phc::R_over_c_p<real_t>(rt_0) / phc::R_d_over_c_pd<real_t>()
	  )) * si::kelvins;

	  arr_z(j) = z / si::metres;
	  arr_rhod(j) = rhod / si::kilograms * si::cubic_metres;
	  arr_zero(j) = 0;
	  arr_rhod_rv(j) = arr_rhod(j) * rt_0;
	  arr_rhod_th(j) = arr_rhod(j) * th / si::kelvins;
	}

	// writing the profiles to the netCDF 
	nvy.putVar(arr_z.data());
	nvrhod.putVar(arr_rhod.data());
	nvrhod_rv.putVar(arr_rhod_rv.data());
	nvrhod_th.putVar(arr_rhod_th.data());

        notice_macro("done.")
      }

      return opts.str();
    } 
  };
};
