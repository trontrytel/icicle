/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#include <libmpdata++/bcond/cyclic_2d.hpp>
#include <libmpdata++/concurr/threads.hpp>

#include <boost/assign/ptr_map_inserter.hpp>  // for 'ptr_map_insert()'

#include "kin_cloud_2d_blk_1m.hpp"
#include "kin_cloud_2d_blk_2m.hpp"
#include "kin_cloud_2d_lgrngn.hpp"

#include "icmw8_case1.hpp" // 8th ICMW case 1 by Wojciech Grabowski)

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
namespace po = boost::program_options;

#include <boost/exception/all.hpp>



// some globals for option handling
int ac;
char** av; // TODO: write it down to a file as in icicle ... write the default (i.e. not specified) values as well!
po::options_description opts_main("General options"); 


void handle_opts(
  po::options_description &opts_micro,
  po::variables_map &vm
)
{
  opts_main.add(opts_micro);
  po::store(po::parse_command_line(ac, av, opts_main), vm); // could be exchanged with a config file parser

  // hendling the "help" option
  if (vm.count("help"))
  {
    std::cout << opts_main;
    exit(EXIT_SUCCESS);
  }
  po::notify(vm); // includes checks for required options
}


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
  po::options_description opts("Single-moment bulk microphysics options"); 
  opts.add_options()
    ("cevp", po::value<bool>()->default_value(true) , "cloud water evaporation (1=on, 0=off)")
    ("revp", po::value<bool>()->default_value(true) , "rain water evaporation (1=on, 0=off)")
    ("conv", po::value<bool>()->default_value(true) , "conversion of cloud water into rain (1=on, 0=off)")
    ("clct", po::value<bool>()->default_value(true) , "cloud water collection by rain (1=on, 0=off)")
    ("sedi", po::value<bool>()->default_value(true) , "rain water sedimentation (1=on, 0=off)")
//TODO: venti
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
  po::options_description opts("Double-moment bulk microphysics options"); 
  opts.add_options()
    ("acti", po::value<bool>()->default_value(true) , "TODO (on/off)")
    ("cond", po::value<bool>()->default_value(true) , "TODO (on/off)")
    ("accr", po::value<bool>()->default_value(true) , "TODO (on/off)")
    ("acnv", po::value<bool>()->default_value(true) , "TODO (on/off)")
    ("turb", po::value<bool>()->default_value(true) , "TODO (on/off)")
    ("sedi", po::value<bool>()->default_value(true) , "TODO (on/off)")
//TODO: venti
  ;
  po::variables_map vm;
  handle_opts(opts, vm);

  // Morrison and Grabowski 2007 scheme options
  params.cloudph_opts.acti = vm["acti"].as<bool>();
  params.cloudph_opts.cond = vm["cond"].as<bool>();
  params.cloudph_opts.accr = vm["accr"].as<bool>();
  params.cloudph_opts.acnv = vm["acnv"].as<bool>();
  params.cloudph_opts.turb = vm["turb"].as<bool>();
  params.cloudph_opts.sedi = vm["sedi"].as<bool>();
}



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
  ;
  po::variables_map vm;
  handle_opts(opts, vm);
      
  thrust_real_t kappa = .5; // TODO!

  //typename libcloudphxx::lgrngn::opts_t<thrust_real_t>::dry_distros_t dry_distros;
/*
  prtcls.reset(libcloudphxx::lgrngn::factory<thrust_real_t>::make(backend,
    vm["sd_conc_mean"].as<thrust_real_t>(), 
    dry_distros,
    nx, params.dx, nz, params.dz
  ));
*/
  std::string backend_str = vm["backend"].as<std::string>();
  if (backend_str == "CUDA") params.backend = libcloudphxx::lgrngn::cuda;
  else if (backend_str == "OpenMP") params.backend = libcloudphxx::lgrngn::omp;
  else if (backend_str == "serial") params.backend = libcloudphxx::lgrngn::cpp;

  params.cloudph_opts.sd_conc_mean = vm["sd_conc_mean"].as<thrust_real_t>();;
  params.cloudph_opts.nx = nx;
  params.cloudph_opts.nz = nz;
  boost::assign::ptr_map_insert<
    icmw8_case1::log_dry_radii<thrust_real_t> // value type
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
}



// model run logic - the same for any microphysics
template <class solver_t>
void run(int nx, int nz, int nt, const std::string &outfile, const int &outfreq)
{
  // instantiation of structure containing simulation parameters
  typename solver_t::params_t p;

  // output and simulation parameters
  p.outfile = outfile;
  p.outfreq = outfreq;
  icmw8_case1::setopts(p, nx, nz);
  setopts_micro<solver_t>(p, nx, nz, nt);

  // solver instantiation
  concurr::threads<solver_t, bcond::cyclic, bcond::cyclic> slv(nx, nz, p);

  // initial condition
  icmw8_case1::intcond(slv);

  // timestepping
  slv.advance(nt);
}



// all starts here with handling general options 
int main(int argc, char** argv)
{
  // making argc and argv global
  ac = argc;
  av = argv;

  const int n_iters = 2; // TODO: where to put such stuff? n_iters should be a param!

  try
  {
    opts_main.add_options()
      ("micro", po::value<std::string>()->required(), "one of: blk_1m, blk_2m, lgrngn")
      ("nx", po::value<int>()->default_value(32) , "grid cell count in horizontal")
      ("nz", po::value<int>()->default_value(32) , "grid cell count in vertical")
      ("nt", po::value<int>()->default_value(500) , "timestep count")
      ("outfile", po::value<std::string>()->required(), "output file name (netCDF-compatible HDF5)")
      ("outfreq", po::value<int>()->required(), "output rate (timestep interval)")
      ("help", "produce a help message (see also --micro X --help)")
    ;
    po::variables_map vm;
    po::store(po::command_line_parser(ac, av).options(opts_main).allow_unregistered().run(), vm); // ignores unknown

    // hendling the "help" option
    if (ac == 1 || (vm.count("help") && !vm.count("micro"))) 
    {
      std::cout << opts_main;
      exit(EXIT_SUCCESS);
    }

    // checking if all required options present
    po::notify(vm); 

    // handling nx, nz, nt options
    int 
      nx = vm["nx"].as<int>(),
      nz = vm["nz"].as<int>(),
      nt = vm["nt"].as<int>();

    // handling the "outfile" option
    std::string outfile = vm["outfile"].as<std::string>();
    int outfreq = vm["outfreq"].as<int>();

    // handling the "micro" option
    std::string micro = vm["micro"].as<std::string>();
    if (micro == "blk_1m")
    {
      struct ix { enum {rhod_th, rhod_rv, rhod_rc, rhod_rr}; };
      run<kin_cloud_2d_blk_1m<icmw8_case1::real_t, n_iters, ix>>(nx, nz, nt, outfile, outfreq);
    }
    else
    if (micro == "blk_2m")
    {
      struct ix { enum {rhod_th, rhod_rv, rhod_rc, rhod_rr, rhod_nc, rhod_nr}; };
      run<kin_cloud_2d_blk_2m<icmw8_case1::real_t, n_iters, ix>>(nx, nz, nt, outfile, outfreq);
    }
    else 
    if (micro == "lgrngn")
    {
      struct ix { enum {rhod_th, rhod_rv}; };
      run<kin_cloud_2d_lgrngn<icmw8_case1::real_t, n_iters, ix>>(nx, nz, nt, outfile, outfreq);
    }
    else BOOST_THROW_EXCEPTION(
      po::validation_error(
        po::validation_error::invalid_option_value, micro, "micro" 
      )
    );
  }
  catch (std::exception &e)
  {
    std::cerr << boost::current_exception_diagnostic_information();
    exit(EXIT_FAILURE);
  }
}
