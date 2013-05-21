/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#pragma once

#  include "opt.hpp"
#  include "grd.hpp"
#  include "ini_func_boxcar.hpp"
#  include "ini_netcdf.hpp"
#  include "ini_func_gauss.hpp"

inline void opt_ini_desc(po::options_description &desc)
{
  desc.add_options()
    ("ini", po::value<string>(), "initical condition: boxcar")

    ("ini.boxcar.ax", po::value<string>()->default_value("0"), "ax [m]")
    ("ini.boxcar.bx", po::value<string>(), "bx [m]")
    ("ini.boxcar.ay", po::value<string>()->default_value("-1"), "ay [m]")
    ("ini.boxcar.by", po::value<string>()->default_value("+1"), "by [m]")
    ("ini.boxcar.az", po::value<string>()->default_value("-1"), "az [m]")
    ("ini.boxcar.bz", po::value<string>()->default_value("+1"), "bz [m]")
    ("ini.boxcar.A", po::value<string>()->default_value("1"), "A [1]")
    ("ini.boxcar.A0", po::value<string>()->default_value("0"), "A0 [1]")

    ("ini.gauss.A", po::value<string>(), "A [1]")
    ("ini.gauss.A0", po::value<string>()->default_value("0"), "A0 [1]")
    ("ini.gauss.x0", po::value<string>(), "x0 [m]")
    ("ini.gauss.y0", po::value<string>()->default_value(".5"), "y0 [m]")
    ("ini.gauss.z0", po::value<string>()->default_value(".5"), "z0 [m]")
    ("ini.gauss.sx", po::value<string>(), "sx [m]")
    ("ini.gauss.sy", po::value<string>()->default_value("1"), "sy [m]")
    ("ini.gauss.sz", po::value<string>()->default_value("1"), "sz [m]")

    ("ini.netcdf.file", po::value<string>()->default_value(""), "input filename")
  ;
}

template <typename real_t>
ini<real_t> *opt_ini(const po::variables_map& vm, const grd<real_t> &grid)
{
  string initype= vm.count("ini") ? vm["ini"].as<string>() : "<unspecified>";
  if (initype == "boxcar")
  {
    quantity<si::length, real_t> 
      ax = real_cast<real_t>(vm, "ini.boxcar.ax") * si::metres,
      bx = real_cast<real_t>(vm, "ini.boxcar.bx") * si::metres,
      ay = real_cast<real_t>(vm, "ini.boxcar.ay") * si::metres,
      by = real_cast<real_t>(vm, "ini.boxcar.by") * si::metres,
      az = real_cast<real_t>(vm, "ini.boxcar.az") * si::metres,
      bz = real_cast<real_t>(vm, "ini.boxcar.bz") * si::metres;
    quantity<si::dimensionless, real_t>
      A  = real_cast<real_t>(vm, "ini.boxcar.A"),
      A0 = real_cast<real_t>(vm, "ini.boxcar.A0");
    return new ini_func_boxcar<real_t>(grid, ax, bx, ay, by, az, bz, A, A0);
  }
#  ifdef USE_NETCDF
  else if (initype == "netcdf")
  {
    return new ini_netcdf<real_t>(grid, vm["ini.netcdf.file" ].as<string>());
  }
#  endif
  else if (initype == "gauss")
  {
    quantity<si::dimensionless, real_t>
      A = real_cast<real_t>(vm, "ini.gauss.A"),
      A0 = real_cast<real_t>(vm, "ini.gauss.A0");
    quantity<si::length, real_t>
      x0 = real_cast<real_t>(vm, "ini.gauss.x0") * si::metres,
      y0 = real_cast<real_t>(vm, "ini.gauss.y0") * si::metres,
      z0 = real_cast<real_t>(vm, "ini.gauss.z0") * si::metres,
      sx = real_cast<real_t>(vm, "ini.gauss.sx") * si::metres,
      sy = real_cast<real_t>(vm, "ini.gauss.sy") * si::metres,
      sz = real_cast<real_t>(vm, "ini.gauss.sz") * si::metres;
    return new ini_func_gauss<real_t>(grid, A, A0, x0, y0, z0, sx, sy, sz);
  }
  else error_macro("unsupported initial condition: " << initype)
}
