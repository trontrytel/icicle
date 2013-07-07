/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#include <libmpdata++/bcond/cyclic_2d.hpp>
#include <libmpdata++/concurr/threads.hpp>
#include <libmpdata++/output/gnuplot.hpp>

#include "kin_cloud_2d_blk_1m.hpp"
#include "kin_cloud_2d_blk_2m.hpp"
#include "kin_cloud_2d_lgrngn.hpp"

#include "icmw8_case1.hpp" // 8th ICMW case 1 by Wojciech Grabowski)

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
namespace po = boost::program_options;

#include <boost/exception/all.hpp>
struct error: virtual boost::exception, virtual std::exception { }; 



// some globals for option handling
int ac;
char** av; // TODO: write it down to a file as in icicle ... write the default (i.e. not specified) values as well!



// simulation and output parameters for micro=blk_1m
template <class solver_t>
void setopts(
  typename solver_t::params_t &params, 
  int nt,
  typename std::enable_if<std::is_same<
    decltype(solver_t::params_t::cloudph_opts),
    libcloudphxx::blk_1m::opts<typename solver_t::real_t>
  >::value>::type* = 0
)
{
  po::options_description opts("Single-moment bulk microphysics options"); 
  opts.add_options()
    ("micro", po::value<std::string>()->required(), "blk_1m")
    ("cevp", po::value<bool>()->default_value(true) , "cloud water evaporation (on/off)")
    ("revp", po::value<bool>()->default_value(true) , "rain water evaporation (on/off)")
    ("conv", po::value<bool>()->default_value(true) , "conversion of cloud water into rain (on/off)")
    ("clct", po::value<bool>()->default_value(true) , "cloud water collection by rain (on/off)")
    ("sedi", po::value<bool>()->default_value(true) , "rain water sedimentation (on/off)")
//TODO: venti
    ("help", "produce a help message")
  ;
  po::variables_map vm;
  po::store(po::parse_command_line(ac, av, opts), vm); // could be exchanged with a config file parser

  // hendling the "help" option
  if (vm.count("help")) 
  {
    std::cout << opts;
    exit(EXIT_SUCCESS);
  }

  using ix = typename solver_t::ix;

  params.outfreq = nt / 50;  // TODO: into general options?
  //params.gnuplot_zrange = p.gnuplot_cbrange = "[.5:2.5]"; // TODO: per variable!
  params.gnuplot_view = "map";
  params.gnuplot_output = "output/figure_%s_%d.svg";
  params.outvars = 
  {
    {ix::rhod_th, {.name = "\\rho_d \\theta", .unit = "kg/m^{-3} K"}},
    {ix::rhod_rv, {.name = "\\rho_v", .unit = "kg/m^{-3}"}},
    {ix::rhod_rc, {.name = "\\rho_c", .unit = "kg/m^{-3}"}},
    {ix::rhod_rr, {.name = "\\rho_r", .unit = "kg/m^{-3}"}}
  };

  // Kessler scheme options
  params.cloudph_opts.cevp = vm["cevp"].as<bool>();
  params.cloudph_opts.revp = vm["revp"].as<bool>();
  params.cloudph_opts.conv = vm["conv"].as<bool>();
  params.cloudph_opts.clct = vm["clct"].as<bool>();
  params.cloudph_opts.sedi = vm["sedi"].as<bool>();
}



// simulation and output parameters for micro=blk_2m
template <class solver_t>
void setopts(
  typename solver_t::params_t &params, 
  int nt,
  typename std::enable_if<std::is_same<
    decltype(solver_t::params_t::cloudph_opts),
    libcloudphxx::blk_2m::opts<typename solver_t::real_t>
  >::value>::type* = 0
)
{
  po::options_description opts("Double-moment bulk microphysics options"); 
  opts.add_options()
    ("micro", po::value<std::string>()->required(), "blk_2m")
    ("acti", po::value<bool>()->default_value(true) , "TODO (on/off)")
    ("cond", po::value<bool>()->default_value(true) , "TODO (on/off)")
    ("accr", po::value<bool>()->default_value(true) , "TODO (on/off)")
    ("acnv", po::value<bool>()->default_value(true) , "TODO (on/off)")
    ("turb", po::value<bool>()->default_value(true) , "TODO (on/off)")
    ("sedi", po::value<bool>()->default_value(true) , "TODO (on/off)")
//TODO: venti
    ("help", "produce a help message")
  ;
  po::variables_map vm;
  po::store(po::parse_command_line(ac, av, opts), vm); // could be exchanged with a config file parser

  // hendling the "help" option
  if (vm.count("help")) 
  {
    std::cout << opts;
    exit(EXIT_SUCCESS);
  }

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
void setopts(
  typename solver_t::params_t &params, 
  int nt,
  typename std::enable_if<std::is_same<
    decltype(solver_t::params_t::cloudph_opts),
    libcloudphxx::lgrngn::opts<typename solver_t::real_t>
  >::value>::type* = 0
)
{
}



// model run logic - the same for any microphysics
template <class solver_t>
void run(int nx, int nz, int nt)
{
  // instantiation of structure containing simulation parameters
  typename solver_t::params_t p;

  // output and simulation parameters
  setopts<solver_t>(p, nt);
  icmw8_case1::setopts(p, nz);

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
    po::options_description opts("General options"); 
    opts.add_options()
      ("micro", po::value<std::string>()->required(), "one of: blk_1m, blk_2m, lgrngn")
      ("nx", po::value<int>()->default_value(32) , "grid cell count in horizontal")
      ("nz", po::value<int>()->default_value(32) , "grid cell count in vertical")
      ("nt", po::value<int>()->default_value(500) , "timestep count")
      ("help", "produce a help message (see also --micro X --help)")
    ;
    po::variables_map vm;
    po::store(po::command_line_parser(ac, av).options(opts).allow_unregistered().run(), vm); // ignores unknown

    // hendling the "help" option
    if (ac == 1 || (vm.count("help") && !vm.count("micro"))) 
    {
      std::cout << opts;
      exit(EXIT_SUCCESS);
    }

    // checking if all required options present
    po::notify(vm); 

    // handling nx, nz, nt options
    int 
      nx = vm["nx"].as<int>(),
      nz = vm["nz"].as<int>(),
      nt = vm["nt"].as<int>();

    // handling the "micro" option
    std::string micro = vm["micro"].as<std::string>();
    if (micro == "blk_1m")
    {
      struct ix { enum {rhod_th, rhod_rv, rhod_rc, rhod_rr}; };
      run<output::gnuplot<kin_cloud_2d_blk_1m<icmw8_case1::real_t, n_iters, solvers::strang, ix>>>(nx, nz, nt);
    }
    else
    if (micro == "blk_2m")
    {
      struct ix { enum {rhod_th, rhod_rv, rhod_rc, rhod_rr, rhod_nc, rhod_nr}; };
      run<output::gnuplot<kin_cloud_2d_blk_2m<icmw8_case1::real_t, n_iters, solvers::strang, ix>>>(nx, nz, nt);
    }
    else 
    if (micro == "lgrngn")
    {
      struct ix { enum {rhod_th, rhod_rv}; };
      run<output::gnuplot<kin_cloud_2d_lgrngn<icmw8_case1::real_t, n_iters, solvers::strang, ix>>>(nx, nz, nt);
    }
    else BOOST_THROW_EXCEPTION(
      po::validation_error(
        po::validation_error::invalid_option_value, micro, "micro" 
      )
    );
  }
  catch (exception &e)
  {
    std::cerr << boost::current_exception_diagnostic_information();
    exit(EXIT_FAILURE);
  }
}
