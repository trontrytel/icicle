/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011 - February 2012
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef OPT_VEL_HPP
#  define OPT_VEL_HPP

#  include "opt.hpp"
#  include "grd.hpp"
#  include "vel_func_uniform.hpp"
#  include "vel_func_rasinski.hpp"
#  include "vel_func_test.hpp"

inline void opt_vel_desc(po::options_description &desc)
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
vel<real_t> *opt_vel(const po::variables_map& vm, grd<real_t> *grid)
{
  string veltype= vm.count("vel") ? vm["vel"].as<string>() : "<unspecified>";
  if (veltype == "uniform")
  {
    quantity<si::velocity, real_t> 
      u = real_cast<real_t>(vm, "vel.uniform.u") * si::metres / si::seconds,
      v = real_cast<real_t>(vm, "vel.uniform.v") * si::metres / si::seconds,
      w = real_cast<real_t>(vm, "vel.uniform.w") * si::metres / si::seconds;
    return new vel_func_uniform<real_t>(u, v, w);
  }
  else if (veltype == "rasinski")
  {
    if (!vm.count("vel.rasinski.Z_clb") || !vm.count("vel.rasinski.Z_top") || !vm.count("vel.rasinski.A"))
      error_macro("vel.rasinski.[Z_clb,Z_top,A] must be specified")
    quantity<si::length, real_t> 
      Z_clb = real_cast<real_t>(vm, "vel.rasinski.Z_clb") * si::metres,
      Z_top = real_cast<real_t>(vm, "vel.rasinski.Z_top") * si::metres,
      X = real_t(grid->nx()) * grid->dx(); // TODO nx+1 dla Arakawa-C ...
    quantity<velocity_times_length, real_t>
      A = real_cast<real_t>(vm, "vel.rasinski.A") * si::metres * si::metres / si::seconds;
    return new vel_func_rasinski<real_t>(X, Z_clb, Z_top, A);
  }
  else if (veltype == "test")
  {
    if (!vm.count("vel.test.omega")) error_macro("vel.test.omega must be specified")
    quantity<si::frequency, real_t> omega = real_cast<real_t>(vm, "vel.test.omega") / si::seconds;
    quantity<si::velocity, real_t> v = real_cast<real_t>(vm, "vel.test.v") * si::metres / si::seconds;
    return new vel_func_test<real_t>(omega, real_t(.5 * grid->nx()) * grid->dx(), real_t(.5 * grid->nz()) * grid->dz(), v);
  }
  else error_macro("unsupported velocity field type: " << veltype)
}

#endif
