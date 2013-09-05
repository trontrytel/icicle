/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include "opts_common.hpp"
#include "kin_cloud_2d_blk_2m.hpp"


// simulation and output parameters for micro=blk_2m
template <class solver_t>
void setopts_micro(
  typename solver_t::params_t &params, 
  int nx, int nz, int nt,
  typename std::enable_if<std::is_same<
    decltype(solver_t::params_t::cloudph_opts),
    libcloudphxx::blk_2m::opts_t<typename solver_t::real_t>
  >::value>::type* = 0
)
{
std::cerr << "setopts_blk_2m" << std::endl;
  po::options_description opts("Double-moment bulk microphysics options"); 
  opts.add_options()
    ("acti", po::value<bool>()->default_value(true) , "TODO (on/off)")
    ("cond", po::value<bool>()->default_value(true) , "TODO (on/off)")
    ("accr", po::value<bool>()->default_value(true) , "TODO (on/off)")
    ("acnv", po::value<bool>()->default_value(true) , "TODO (on/off)")
    ("sedi", po::value<bool>()->default_value(true) , "TODO (on/off)")
    ("mean_rd", po::value<setup::real_t>()->default_value(setup::mean_rd / si::metres) , "TODO")
    ("sdev_rd", po::value<setup::real_t>()->default_value(setup::sdev_rd) , "TODO")
    ("N_stp",   po::value<setup::real_t>()->default_value(setup::N_stp * si::cubic_metres  ) , "TODO")
    ("chem_b",  po::value<setup::real_t>()->default_value(setup::chem_b ) , "TODO")
//TODO: venti
  ;
  po::variables_map vm;
  handle_opts(opts, vm);

  // Morrison and Grabowski 2007 scheme options
  params.cloudph_opts.acti = vm["acti"].as<bool>();
  params.cloudph_opts.cond = vm["cond"].as<bool>();
  params.cloudph_opts.accr = vm["accr"].as<bool>();
  params.cloudph_opts.acnv = vm["acnv"].as<bool>();
  params.cloudph_opts.sedi = vm["sedi"].as<bool>();
  params.cloudph_opts.mean_rd = vm["mean_rd"].as<setup::real_t>();
  params.cloudph_opts.sdev_rd = vm["sdev_rd"].as<setup::real_t>();
  params.cloudph_opts.N_stp   = vm["N_stp"].as<setup::real_t>();
  params.cloudph_opts.chem_b  = vm["chem_b"].as<setup::real_t>();

  // output variables
  params.outvars = {
    // <TODO>: make it common among all three micro?
    {solver_t::ix::rhod_th, {"rhod_th", "[K kg m-3]"}},
    {solver_t::ix::rhod_rv, {"rhod_rv", "[kg m-3]"}},
    // </TODO>
    {solver_t::ix::rhod_rc, {"rhod_rc", "[kg m-3]"}},
    {solver_t::ix::rhod_rr, {"rhod_rr", "[kg m-3]"}},
    {solver_t::ix::rhod_nc, {"rhod_nc", "[m-3]"}},
    {solver_t::ix::rhod_nr, {"rhod_nr", "[m-3]"}}
  };
}
