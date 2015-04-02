/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#include <libmpdata++/bcond/cyclic_2d.hpp>
#include <libmpdata++/bcond/open_2d.hpp>
#include <libmpdata++/concurr/boost_thread.hpp> // not to conflict with OpenMP used via Thrust in libcloudph++

#include "icmw8_case1.hpp" // 8th ICMW case 1 by Wojciech Grabowski)
namespace setup = icmw8_case1;

#include "opts_blk_1m.hpp"
#include "opts_blk_2m.hpp"
#include "opts_lgrngn.hpp"

// exception handling
#include <boost/exception/all.hpp>

#include "panic.hpp"

// model run logic - the same for any microphysics
template <class solver_t>
void run(int nx, int nz, int nt, const std::string &outdir, const int &outfreq, int spinup)
{
  // instantiation of structure containing simulation parameters
  typename solver_t::rt_params_t p;

  // output and simulation parameters
  p.grid_size = {nx, nz};
  p.outdir = outdir;
  p.outfreq = outfreq;
  p.spinup = spinup;
  setup::setopts(p, nx, nz);
  setopts_micro<solver_t>(p, nx, nz, nt);

  // solver instantiation
  concurr::boost_thread<solver_t, 
    bcond::cyclic, bcond::cyclic,
    bcond::open,   bcond::open 
  > slv(p);

  // initial condition
  setup::intcond(slv);

  // setup panic pointer and the signal handler
  panic = slv.panic_ptr();
  set_sigaction();
 
  // timestepping
  slv.advance(nt);
}


// libmpdata++'s compile-time parameters
struct ct_params_common : ct_params_default_t
{
  using real_t = setup::real_t;
  enum { n_dims = 2 };
  enum { opts = opts::nug | opts::fct }; 
  enum { rhs_scheme = solvers::euler_b };
};


// all starts here with handling general options 
int main(int argc, char** argv)
{
  // making argc and argv global
  ac = argc;
  av = argv;

  try
  {
    // note: all options should have default values here to make "--micro=? --help" work
    opts_main.add_options()
      ("micro", po::value<std::string>()->required(), "one of: blk_1m, blk_2m, lgrngn")
      ("nx", po::value<int>()->default_value(101) , "grid cell count in horizontal")
      ("nz", po::value<int>()->default_value(51) , "grid cell count in vertical")
      ("nt", po::value<int>()->default_value(3600) , "timestep count") // TODO!
      ("outdir", po::value<std::string>(), "output file name (netCDF-compatible HDF5)")
      ("outfreq", po::value<int>(), "output rate (timestep interval)")
      ("spinup", po::value<int>()->default_value(2400) , "number of initial timesteps during which rain formation is to be turned off") //TODO!
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
    
    // handling outdir && outfreq
    std::string outdir; 
    int outfreq;
    if (!vm.count("help"))
    {
      if (!vm.count("outdir")) BOOST_THROW_EXCEPTION(po::required_option("outdir"));
      if (!vm.count("outfreq")) BOOST_THROW_EXCEPTION(po::required_option("outfreq"));
      outdir = vm["outdir"].as<std::string>();
      outfreq = vm["outfreq"].as<int>();
    }

    // handling nx, nz, nt options
    int 
      nx = vm["nx"].as<int>(),
      nz = vm["nz"].as<int>(),
      nt = vm["nt"].as<int>(),
      spinup = vm["spinup"].as<int>();

    // handling the "micro" option
    std::string micro = vm["micro"].as<std::string>();
    if (micro == "blk_1m")
    {
      // libmpdata++'s compile-time parameters
      struct ct_params_t : ct_params_common
      {
	enum { n_eqns = 4 };
        struct ix { enum {th, rv, rc, rr}; };
        enum { hint_norhs = opts::bit(ix::th) | opts::bit(ix::rv) }; // only through adjustments
      };
      run<kin_cloud_2d_blk_1m<ct_params_t>>(nx, nz, nt, outdir, outfreq, spinup);
    }
    else
    if (micro == "blk_2m")
    {
      struct ct_params_t : ct_params_common
      {
	enum { n_eqns = 6 };
	struct ix { enum {th, rv, rc, rr, nc, nr}; }; 

        static constexpr int hint_scale(const int &e) 
	{
	  return 
            e == ix::nc ?  24 : // 1.7e7
            e == ix::nr ?  17 : // 1.3e5 
            e == ix::rc ? -14 : // 1.6e4
            e == ix::rr ? -14 : // 1.6e4
            0;
	}
      };
      run<kin_cloud_2d_blk_2m<ct_params_t>>(nx, nz, nt, outdir, outfreq, spinup);
    }
    else 
    if (micro == "lgrngn")
    {
      struct ct_params_t : ct_params_common
      {
	enum { n_eqns = 2 };
	struct ix { enum {th, rv}; };
        enum { hint_norhs = opts::bit(ix::th) | opts::bit(ix::rv) }; // only through adjustments
      };
      run<kin_cloud_2d_lgrngn<ct_params_t>>(nx, nz, nt, outdir, outfreq, spinup);
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
