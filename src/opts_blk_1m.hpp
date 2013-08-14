/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include "opts_common.hpp"
#include "kin_cloud_2d_blk_1m.hpp"

// simulation and output parameters for micro=blk_1m
template <class solver_t>
void setopts_micro(
  typename solver_t::params_t &params, 
  int nx, int nz, int nt,
  typename std::enable_if<std::is_same<
    decltype(solver_t::params_t::cloudph_opts),
    libcloudphxx::blk_1m::opts_t<typename solver_t::real_t>
  >::value>::type* = 0
)
{
std::cerr << "setopts_blk_1m" << std::endl;
  po::options_description opts("Single-moment bulk microphysics options"); 
  opts.add_options()
    ("cevp", po::value<bool>()->default_value(true) , "cloud water evaporation (1=on, 0=off)")
    ("revp", po::value<bool>()->default_value(true) , "rain water evaporation (1=on, 0=off)")
    ("conv", po::value<bool>()->default_value(true) , "autoconversion of cloud water into rain (1=on, 0=off)")
    ("clct", po::value<bool>()->default_value(true) , "cloud water collection by rain (1=on, 0=off)")
    ("sedi", po::value<bool>()->default_value(true) , "rain water sedimentation (1=on, 0=off)")
//TODO: venti, autoconv_threshold
  ;
  po::variables_map vm;
  handle_opts(opts, vm);

  // Kessler scheme options
  params.cloudph_opts.cevp = vm["cevp"].as<bool>();
  params.cloudph_opts.revp = vm["revp"].as<bool>();
  params.cloudph_opts.conv = vm["conv"].as<bool>();
  params.cloudph_opts.clct = vm["clct"].as<bool>();
  params.cloudph_opts.sedi = vm["sedi"].as<bool>();

  // output variables
  params.outvars = {
    // <TODO>: make it common among all three micro?
    {solver_t::ix::rhod_th, {"rhod_th", "[K kg m-3]"}},
    {solver_t::ix::rhod_rv, {"rhod_rv", "[kg m-3]"}},
    // </TODO>
    {solver_t::ix::rhod_rc, {"rhod_rc", "[kg m-3]"}},
    {solver_t::ix::rhod_rr, {"rhod_rr", "[kg m-3]"}}
  };
}
