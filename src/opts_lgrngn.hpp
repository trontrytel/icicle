/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#include <boost/assign/ptr_map_inserter.hpp>  // for 'ptr_map_insert()'

#include "opts_common.hpp"
#include "kin_cloud_2d_lgrngn.hpp"

// string parsing
#include <boost/spirit/include/qi.hpp>    
#include <boost/fusion/adapted/std_pair.hpp> 

// simulation and output parameters for micro=lgrngn
template <class solver_t>
void setopts_micro(
  typename solver_t::params_t &params, 
  int nx, int nz, int nt,
  typename std::enable_if<std::is_same<
    decltype(solver_t::params_t::cloudph_opts),
    libcloudphxx::lgrngn::opts_t<typename solver_t::real_t>
  >::value>::type* = 0
)
{
  using thrust_real_t = double; //float; // TODO: option, warning, ...?  (if nvcc downgraded real_t=double to float)

  po::options_description opts("Lagrangian microphysics options"); 
  opts.add_options()
    ("backend", po::value<std::string>()->required() , "backend (one of: CUDA, OpenMP, serial)")
    ("sd_conc_mean", po::value<thrust_real_t>()->required() , "mean super-droplet concentration per grid cell (int)")
    // processes
    ("adve", po::value<bool>()->default_value(true) , "particle advection             (1=on, 0=off)")
    ("sedi", po::value<bool>()->default_value(true) , "particle sedimentation         (1=on, 0=off)")
    ("rcyc", po::value<bool>()->default_value(true) , "particle recycling             (1=on, 0=off)")
    ("cond", po::value<bool>()->default_value(true) , "particle condensational growth (1=on, 0=off)")
    ("coal", po::value<bool>()->default_value(true) , "particle collisional growth    (1=on, 0=off)")
    // 
    ("out_md0", po::value<bool>()->default_value(""),            "dry radius ranges for 0-th moment output (r1:r2;r2:r3;...)")
    ("out_mw0", po::value<bool>()->default_value(".5e-6:25e-6"), "wet      --||--       0-th     --||--")
    ("out_mw1", po::value<bool>()->default_value(".5e-6:25e-6"), "wet      --||--       1-st     --||--")
    ("out_mw2", po::value<bool>()->default_value(".5e-6:25e-6"), "wet      --||--       2-nd     --||--")
    ("out_mw3", po::value<bool>()->default_value(".5e-6:25e-6"), "wet      --||--       3-rd     --||--")
    ("out_mw6", po::value<bool>()->default_value(""),            "wet      --||--       6-th     --||--")
  ;
  po::variables_map vm;
  handle_opts(opts, vm);
      
  thrust_real_t kappa = .5; // TODO!!!

  std::string backend_str = vm["backend"].as<std::string>();
  if (backend_str == "CUDA") params.backend = libcloudphxx::lgrngn::cuda;
  else if (backend_str == "OpenMP") params.backend = libcloudphxx::lgrngn::omp;
  else if (backend_str == "serial") params.backend = libcloudphxx::lgrngn::cpp;

  params.cloudph_opts.sd_conc_mean = vm["sd_conc_mean"].as<thrust_real_t>();;
  params.cloudph_opts.nx = nx;
  params.cloudph_opts.nz = nz;
  boost::assign::ptr_map_insert<
    setup::log_dry_radii<thrust_real_t> // value type
  >(
    params.cloudph_opts.dry_distros // map
  )(
    kappa // key
  );

  // output variables
  params.outvars = {
    // <TODO>: make it common among all three micro?
    {solver_t::ix::rhod_th, {"rhod_th", "[K kg m-3]"}},
    {solver_t::ix::rhod_rv, {"rhod_rv", "[kg m-3]"}}
    // </TODO>
  };

  // process toggling
  params.cloudph_opts.adve = vm["adve"].as<bool>();
  params.cloudph_opts.sedi = vm["sedi"].as<bool>();
  params.cloudph_opts.rcyc = vm["rcyc"].as<bool>();
  params.cloudph_opts.cond = vm["cond"].as<bool>();
  params.cloudph_opts.coal = vm["coal"].as<bool>();

  // output moments
  
}
