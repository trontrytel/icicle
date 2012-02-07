/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef OPT_INI_HPP
#  define OPT_INI_HPP

#  include "opt.hpp"
#  include "grd.hpp"
#  include "ini_func_boxcar.hpp"
#  include "ini_func_pks_wwg_1989.hpp"
#  include "ini_func_cone.hpp"
#  include "ini_netcdf.hpp"
#  include "ini_func_gauss.hpp"

inline void opt_ini_desc(po::options_description &desc)
{
  desc.add_options()
    ("ini", po::value<string>(), "initical condition: boxcar, cone, pks_wwg_1989")

    ("ini.boxcar.ax", po::value<string>()->default_value("0"), "ax [m]")
    ("ini.boxcar.bx", po::value<string>(), "bx [m]")
    ("ini.boxcar.ay", po::value<string>()->default_value("-1"), "ay [m]")
    ("ini.boxcar.by", po::value<string>()->default_value("+1"), "by [m]")
    ("ini.boxcar.az", po::value<string>()->default_value("-1"), "az [m]")
    ("ini.boxcar.bz", po::value<string>()->default_value("+1"), "bz [m]")
    ("ini.boxcar.A", po::value<string>()->default_value("1"), "A [1]")
    ("ini.boxcar.A0", po::value<string>()->default_value("0"), "A0 [1]")

    ("ini.cone.h", po::value<string>()->default_value("3.87"), "h [m]")
    ("ini.cone.x0", po::value<string>()->default_value("75"), "x0 [m]")
    ("ini.cone.z0", po::value<string>()->default_value("50"), "z0 [m]")
    ("ini.cone.r", po::value<string>()->default_value("15"), "r [m]")
    ("ini.cone.h0", po::value<string>()->default_value("0"), "h0 [m]")

    ("ini.gauss.A", po::value<string>(), "A [1]")
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
  if (initype == "pks_wwg_1989")
    return new ini_func_pks_wwg_1989<real_t>(grid);
  else if (initype == "boxcar")
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
  else if (initype == "cone")
  {
    quantity<si::length, real_t> 
      h  = real_cast<real_t>(vm, "ini.cone.h" ) * si::metres,
      x0 = real_cast<real_t>(vm, "ini.cone.x0") * si::metres,
      z0 = real_cast<real_t>(vm, "ini.cone.z0") * si::metres,
      r  = real_cast<real_t>(vm, "ini.cone.r" ) * si::metres,
      h0 = real_cast<real_t>(vm, "ini.cone.h0") * si::metres;
    return new ini_func_cone<real_t>(grid, h, x0, z0, r, h0);
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
      A = real_cast<real_t>(vm, "ini.gauss.A");
    quantity<si::length, real_t>
      x0 = real_cast<real_t>(vm, "ini.gauss.x0") * si::metres,
      y0 = real_cast<real_t>(vm, "ini.gauss.y0") * si::metres,
      z0 = real_cast<real_t>(vm, "ini.gauss.z0") * si::metres,
      sx = real_cast<real_t>(vm, "ini.gauss.sx") * si::metres,
      sy = real_cast<real_t>(vm, "ini.gauss.sy") * si::metres,
      sz = real_cast<real_t>(vm, "ini.gauss.sz") * si::metres;
    return new ini_func_gauss<real_t>(grid, A, x0, y0, z0, sx, sy, sz);
  }
  else error_macro("unsupported initial condition: " << initype)
}

#endif
