/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011 - February 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef OPT_VEL_HPP
#  define OPT_VEL_HPP

#  include "opt.hpp"
#  include "grd.hpp"
#  include "vel_func_uniform.hpp"
#  include "vel_func_stream_rasinski.hpp"
#  include "vel_func_test.hpp"
#  include "vel_momeq_extrapol.hpp"

inline void opt_vel_desc(po::options_description &desc)
{
  desc.add_options()
    ("vel", po::value<string>(), "velocity field: uniform, rasinski, test, momeq_extrapol")

    ("vel.uniform.u", po::value<string>()->default_value("0"), "velocity (X) [m/s]")
    ("vel.uniform.v", po::value<string>()->default_value("0"), "velocity (Y) [m/s]")
    ("vel.uniform.w", po::value<string>()->default_value("0"), "velocity (Z) [m/s]")

    ("vel.rasinski.file", po::value<string>(), "netCDF filename (rho(z) profile)")
    ("vel.rasinski.A", po::value<string>(), "amplitude [kg/m/s]")

    ("vel.test.omega", po::value<string>(), "frequency [1/s]")
    ("vel.test.v", po::value<string>()->default_value("0"), "Y velocity [m/s]");
}

template <typename real_t>
vel<real_t> *opt_vel(const po::variables_map& vm, const grd<real_t> &grid)
{
  string veltype= vm.count("vel") ? vm["vel"].as<string>() : "<unspecified>";
  if (veltype == "uniform")
  {
    quantity<si::velocity, real_t> 
      u = real_cast<real_t>(vm, "vel.uniform.u") * si::metres / si::seconds,
      v = real_cast<real_t>(vm, "vel.uniform.v") * si::metres / si::seconds,
      w = real_cast<real_t>(vm, "vel.uniform.w") * si::metres / si::seconds;
    return new vel_func_uniform<real_t>(grid, u, v, w);
  }
  else if (veltype == "test")
  {
    if (!vm.count("vel.test.omega")) error_macro("vel.test.omega must be specified")
    quantity<si::frequency, real_t> omega = real_cast<real_t>(vm, "vel.test.omega") / si::seconds;
    quantity<si::velocity, real_t> v = real_cast<real_t>(vm, "vel.test.v") * si::metres / si::seconds;
    return new vel_func_test<real_t>(grid, omega, real_t(.5 * grid.nx()) * grid.dx(), real_t(.5 * grid.nz()) * grid.dz(), v); // TODO: we pass grid, so part of it could be calculated in the ctor
  }
  else if (veltype == "momeq_extrapol")
  {
    return new vel_momeq_extrapol<real_t>(grid);
  }
  else if (veltype == "rasinski")
  {
    if (!vm.count("vel.rasinski.A") || !vm.count("vel.rasinski.file"))
      error_macro("vel.rasinski.[A,file] must be specified")
    return new vel_func_stream_rasinski<real_t>(grid, vm["vel.rasinski.file" ].as<string>(), 
      real_cast<real_t>(vm, "vel.rasinski.A") * si::kilograms / si::metres / si::seconds
    );
  }
  else error_macro("unsupported velocity field type: " << veltype)
}

#endif
