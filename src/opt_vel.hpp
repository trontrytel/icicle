/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef OPT_VEL_HPP
#  define OPT_VEL_HPP

#  include "opt.hpp"
#  include "grd.hpp"
#  include "vel_uniform.hpp"
#  include "vel_rasinski.hpp"
#  include "vel_test.hpp"

void opt_vel_desc(po::options_description &desc)
{
  desc.add_options()
    ("vel", po::value<string>(), "velocity field: uniform, rasinski")
    ("vel.uniform.u", po::value<string>()->default_value("0"), "velocity (X) [m/s]")
    ("vel.uniform.v", po::value<string>()->default_value("0"), "velocity (Y) [m/s]")
    ("vel.uniform.w", po::value<string>()->default_value("0"), "velocity (Z) [m/s]")
    ("vel.rasinski.Z_clb", po::value<string>(), "cloud base height [m]")
    ("vel.rasinski.Z_top", po::value<string>(), "cloud top height [m]")
    ("vel.rasinski.A", po::value<string>(), "amplitude [m2/s]")
    ("vel.test.omega", po::value<string>(), "frequency [1/s]")
    ("vel.test.v", po::value<string>()->default_value("0"), "Y velocity [m/s]");
}

template <typename real_t>
vel<real_t> *opt_vel(const po::variables_map& vm,
  grd<real_t> *grid, int nx, int ny, int nz)
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
    return new vel_uniform<real_t>(u, v, w);
  }
  else if (veltype == "rasinski")
  {
    if (!vm.count("vel.rasinski.Z_clb") || !vm.count("vel.rasinski.Z_top") || !vm.count("vel.rasinski.A"))
      error_macro("vel.rasinski.[Z_clb,Z_top,A] must be specified")
    quantity<si::length, real_t> 
      Z_clb = boost::lexical_cast<real_t>(vm["vel.rasinski.Z_clb"].as<string>()) * si::metres,
      Z_top = boost::lexical_cast<real_t>(vm["vel.rasinski.Z_top"].as<string>()) * si::metres,
      X = real_t(nx) * grid->dx(); // TODO nx+1 dla Arakawa-C ...
    quantity<velocity_times_length, real_t>
      A = boost::lexical_cast<real_t>(vm["vel.rasinski.A"].as<string>()) * si::metres * si::metres / si::seconds;
    return new vel_rasinski<real_t>(X, Z_clb, Z_top, A);
  }
  else if (veltype == "test")
  {
    if (!vm.count("vel.test.omega")) error_macro("vel.test.omega must be specified")
    quantity<si::frequency, real_t> omega = boost::lexical_cast<real_t>(vm["vel.test.omega"].as<string>()) / si::seconds;
    quantity<si::velocity, real_t> v = boost::lexical_cast<real_t>(vm["vel.test.v"].as<string>()) * si::metres / si::seconds;
    return new vel_test<real_t>(omega, real_t(.5 * nx) * grid->dx(), real_t(.5 * nz) * grid->dz(), v);
  }
  else error_macro("unsupported velocity field type: " << veltype)
}

#endif
