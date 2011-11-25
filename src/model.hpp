/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef MODEL_HPP
#  define MODEL_HPP

#  include "opt.hpp"

#  include "adv_mpdata.hpp"
#  include "adv_leapfrog.hpp"
#  include "adv_lax-wendroff.hpp"

#  include "dom_serial.hpp"
#  include "dom_parallel_openmp.hpp"
#  include "dom_parallel_threads.hpp"

#  include "out_gnuplot.hpp"
#  include "out_netcdf.hpp"

#  include "vel_uniform.hpp"
#  include "vel_2d-xz_rasinski.hpp"
//#  include "vel_2d-xz_test.hpp"

#  include <boost/lexical_cast.hpp>

template <typename real_t>
void model(const po::variables_map& vm) 
{
  // some key parameters 
  if (
    !vm.count("nx") || !vm.count("ny") || !vm.count("nz") || 
    !vm.count("grd.dx") || !vm.count("grd.dy") || !vm.count("grd.dz") ||
    !vm.count("nt") || !vm.count("dt")
  )
    error_macro("nx, ny, nz, grd.dx, grd.dy, grd.dz, nt and dt options are mandatory")
  quantity<si::length, real_t> 
    dx = boost::lexical_cast<real_t>(vm["grd.dx"].as<string>()) * si::metres,
    dy = boost::lexical_cast<real_t>(vm["grd.dy"].as<string>()) * si::metres,
    dz = boost::lexical_cast<real_t>(vm["grd.dz"].as<string>()) * si::metres;
  assert(dx == dy && dy == dz); // TODO
  quantity<si::time, real_t> 
    dt = boost::lexical_cast<real_t>(vm["dt"].as<string>()) * si::seconds;
  int 
    nx = vm["nx"].as<int>(),
    ny = vm["ny"].as<int>(),
    nz = vm["nz"].as<int>();
  unsigned long
    nt = vm["nt"].as<unsigned long>();

  // sanity checks
  if (nx <= 0 || ny <= 0 || nz <= 0 ||
    dx / si::metres <= 0 || dy / si::metres <= 0 || dz / si::metres <= 0 
  ) error_macro("dx, dy, dz, nx, ny, nz must all be >= 0") 

  // grid (FIXME: no deallocation yet)
  grd<real_t> *grid;
  {
    string grdtype = vm.count("grd") ? vm["grd"].as<string>() : "<unspecified>";
    if (grdtype == "arakawa-c")
    {
      grid = new grd_2d_xz_arakawa_c<real_t>(dx, dz);
    }
    else error_macro("unsupported grid type: " << grdtype)
  }

  // advection sheme "allocation" (FIXME: no deallocation yet)
  adv<si::dimensionless, real_t> *advsch, *fllbck = NULL;
  {
    string advscheme = vm.count("adv") ? vm["adv"].as<string>() : "<unspecified>";
    if (advscheme == "leapfrog")
    {
      advsch = new adv_leapfrog<si::dimensionless, real_t>(grid);
      fllbck = new adv_mpdata<si::dimensionless, real_t>(grid, 1);
    }
    else if (advscheme == "mpdata")
    {
      int iord = vm.count("adv.mpdata.iord") ? vm["adv.mpdata.iord"].as<int>() : 2;
      advsch = new adv_mpdata<si::dimensionless, real_t>(grid, iord);
    }
    else error_macro("unsupported advection scheme: " << advscheme)
  }

  // output set-up (FIXME: no deallocation yet)
  out<si::dimensionless, real_t> *output;
  {
    string outtype = vm.count("out") ? vm["out"].as<string>() : "<unspecified>";
    if (outtype == "gnuplot")
      output = new out_gnuplot<si::dimensionless, real_t>();
    else
#  ifdef USE_NETCDF
    if (outtype == "netcdf")
    {
      if (!vm.count("out.netcdf.file")) error_macro("output filename not specified (--out.netcdf.file option)")
      output = new out_netcdf<si::dimensionless, real_t>(vm["out.netcdf.file"].as<string>(), grid, nx, ny, nz);
    }
    else
#  endif
    error_macro("unsupported output type: " << outtype)
  }

  // velocity field // TODO: no deallocation yet!
  vel<real_t> *velocity;
  { 
    string veltype= vm.count("vel") ? vm["vel"].as<string>() : "<unspecified>";
    if (veltype == "uniform")
    {
      if (!vm.count("vel.uniform.u") || !vm.count("vel.uniform.v") || !vm.count("vel.uniform.w"))
        error_macro("vel.uniform.[u,v,w] must be specified")
      quantity<si::velocity, real_t> 
        u = boost::lexical_cast<real_t>(vm["vel.uniform.u"].as<string>()) * si::metres / si::seconds,
        v = boost::lexical_cast<real_t>(vm["vel.uniform.v"].as<string>()) * si::metres / si::seconds,
        w = boost::lexical_cast<real_t>(vm["vel.uniform.w"].as<string>()) * si::metres / si::seconds;
      velocity = new vel_uniform<real_t>(u, v, w);
    }
    else if (veltype == "rasinski")
    {
      if (!vm.count("vel.rasinski.Z_clb") || !vm.count("vel.rasinski.Z_top") || !vm.count("vel.rasinski.A"))
        error_macro("vel.rasinski.[Z_clb,Z_top,A] must be specified")
      quantity<si::length, real_t> 
        Z_clb = boost::lexical_cast<real_t>(vm["vel.rasinski.Z_clb"].as<string>()) * si::metres,
        Z_top = boost::lexical_cast<real_t>(vm["vel.rasinski.Z_top"].as<string>()) * si::metres, 
        X = real_t(nx) * dx; // TODO nx+1 dla Arakawa-C ...
      quantity<velocity_times_length, real_t> 
        A = boost::lexical_cast<real_t>(vm["vel.rasinski.A"].as<string>()) * si::metres * si::metres / si::seconds;
      velocity = new vel_2d_xz_rasinski<real_t>(X, Z_clb, Z_top, A);
    }
//    else if (veltype == "test")
//   {
//      if (!vm.count("vel.test.omega")) error_macro("vel.test.omega must be specified")
//      quantity<si::frequency, real_t> omega = boost::lexical_cast<real_t>(vm["vel.test.omega"].as<string>()) / si::seconds;
//      velocity = new vel_2d_xz_test<real_t>(omega, .5 * nx * dx, .5 * nz * dz);
//    }
    else error_macro("unsupported velocity field type: " << veltype)
  }

  // time integration loop
  // domain set-up (FIXME: no deallocation yet)
  // TODO: rename domain -> solver?
  dom<si::dimensionless, real_t> *domain;
  {
    string domtype = vm.count("dom") ? vm["dom"].as<string>() : "<unspecified>";
    if (domtype == "serial")
      domain = new dom_serial<si::dimensionless, real_t>(fllbck, advsch, output, velocity,
        0, nx - 1, nx, 
        0, ny - 1, ny, 
        0, nz - 1, nz, 
        grid, dt
      );
    else 
    {
      if (!vm.count("nsd")) error_macro("subdomain count not specified (--nsd option)")
      int nsd = vm["nsd"].as<int>();
#  ifdef _OPENMP
      if (domtype == "openmp")
        domain = new dom_parallel_openmp<si::dimensionless, real_t>(
          fllbck, advsch, output, velocity, nx, ny, nz, grid, dt, nsd);
      else 
#  endif
#  ifdef USE_BOOST_THREAD
      if (domtype == "threads")
        domain = new dom_parallel_threads<si::dimensionless, real_t>(
          fllbck, advsch, output, velocity, nx, ny, nz, grid, dt, nsd);
      else 
#  endif
      error_macro("unsupported domain type: " << domtype)
    }
  }

  domain->integ_loop(nt, dt);
}
#endif
