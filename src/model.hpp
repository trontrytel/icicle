/** @file
 *  @author Sylwester Arabas (University of Warsaw) <slayoo@igf.fuw.edu.pl>
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef MODEL_HPP
#  define MODEL_HPP

#  include "opt.hpp"

#  include "adv_upstream.hpp"
#  include "adv_mpdata.hpp"
#  include "adv_leapfrog.hpp"

#  include "dom_serial.hpp"
#  include "dom_parallel_openmp.hpp"
#  include "dom_parallel_threads.hpp"

#  include "out_gnuplot.hpp"
#  include "out_netcdf.hpp"

#  include <boost/lexical_cast.hpp>

template <typename real_t>
void model(const po::variables_map& vm) 
{
  // some key parameters 
  if (
    !vm.count("nx") || !vm.count("ny") || !vm.count("nz") || 
    !vm.count("dx") || !vm.count("dy") || !vm.count("dz") ||
    !vm.count("nt") || !vm.count("dt") || 
    !vm.count("u") || !vm.count("v") || !vm.count("w")
  )
    error_macro("nx, ny, nz, dx, dy, dz, nt, dt and u options are mandatory")
  quantity<si::length, real_t> 
    dx = boost::lexical_cast<real_t>(vm["dx"].as<string>()) * si::metres,
    dy = boost::lexical_cast<real_t>(vm["dy"].as<string>()) * si::metres,
    dz = boost::lexical_cast<real_t>(vm["dz"].as<string>()) * si::metres;
  quantity<si::time, real_t> 
    dt = boost::lexical_cast<real_t>(vm["dt"].as<string>()) * si::seconds;
  quantity<si::velocity, real_t> 
    u = boost::lexical_cast<real_t>(vm["u"].as<string>()) * si::metres / si::seconds,
    v = boost::lexical_cast<real_t>(vm["v"].as<string>()) * si::metres / si::seconds,
    w = boost::lexical_cast<real_t>(vm["w"].as<string>()) * si::metres / si::seconds;
  quantity<si::dimensionless, real_t> 
    Cx = u * dt / dx,
    Cy = v * dt / dx,
    Cz = w * dt / dx;
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

  // advection sheme "allocation" (FIXME: no deallocation yet)
  adv<si::dimensionless, real_t> *advsch, *fllbck = NULL;
  {
    string advscheme = vm.count("adv") ? vm["adv"].as<string>() : "<unspecified>";
    if (advscheme == "upstream") 
      advsch = new adv_upstream<si::dimensionless, real_t>();
    else if (advscheme == "leapfrog")
    {
      advsch = new adv_leapfrog<si::dimensionless, real_t>();
      fllbck = new adv_upstream<si::dimensionless, real_t>();
    }
    else if (advscheme == "mpdata")
    {
      int iord = vm.count("adv.mpdata.iord") ? vm["adv.mpdata.iord"].as<int>() : 2;
      advsch = new adv_mpdata<si::dimensionless, real_t>(iord);
    }
    else error_macro("unsupported advection scheme: " << advscheme);
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
      output = new out_netcdf<si::dimensionless, real_t>(vm["out.netcdf.file"].as<string>(), nx, ny, nz);
    }
    else
#  endif
    error_macro("unsupported output type: " << outtype);
  }

  // domain set-up (FIXME: no deallocation yet)
  dom<si::dimensionless, real_t> *domain;
  {
    string domtype = vm.count("dom") ? vm["dom"].as<string>() : "<unspecified>";
    if (domtype == "serial")
      domain = new dom_serial<si::dimensionless, real_t>(fllbck, advsch, output, 
        0, nx - 1, nx,
        0, ny - 1, ny,
        0, nz - 1, nz
      );
    else 
    {
      if (!vm.count("nsd")) error_macro("subdomain count not specified (--nsd option)");
      int nsd = vm["nsd"].as<int>();
#  ifdef _OPENMP
      if (domtype == "openmp")
        domain = new dom_parallel_openmp<si::dimensionless, real_t>(fllbck, advsch, output, nx, ny, nz, nsd);
      else 
#  endif
#  ifdef USE_BOOST_THREAD
      if (domtype == "threads")
        domain = new dom_parallel_threads<si::dimensionless, real_t>(fllbck, advsch, output, nx, ny, nz, nsd);
      else 
#  endif
      error_macro("unsupported domain type: " << domtype);
    }
  }

  // time integration loop
  domain->integ_loop(nt, Cx, Cy, Cz);
}
#endif
