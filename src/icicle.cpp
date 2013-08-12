/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#include <libmpdata++/bcond/cyclic_2d.hpp>
#include <libmpdata++/concurr/threads.hpp>

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
void run(int nx, int nz, int nt, const std::string &outfile, const int &outfreq)
{
  // instantiation of structure containing simulation parameters
  typename solver_t::params_t p;

  // output and simulation parameters
  p.outfile = outfile;
  p.outfreq = outfreq;
  setup::setopts(p, nx, nz);
  setopts_micro<solver_t>(p, nx, nz, nt);

  // solver instantiation
  concurr::threads<solver_t, bcond::cyclic, bcond::cyclic> slv(nx, nz, p);

  // initial condition
  setup::intcond(slv);

  // setup panic pointer and the signal handler
  panic = slv.panic_ptr();
  set_sigaction();
 
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
    // note: all options should have default values here to make "--micro=? --help" work
    opts_main.add_options()
      ("micro", po::value<std::string>()->required(), "one of: blk_1m, blk_2m, lgrngn")
      ("nx", po::value<int>()->default_value(32) , "grid cell count in horizontal")
      ("nz", po::value<int>()->default_value(32) , "grid cell count in vertical")
      ("nt", po::value<int>()->default_value(500) , "timestep count")
      ("outfile", po::value<std::string>(), "output file name (netCDF-compatible HDF5)")
      ("outfreq", po::value<int>(), "output rate (timestep interval)")
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
    
    // handling outfile && outfreq
    std::string outfile; 
    int outfreq;
    if (!vm.count("help"))
    {
      if (!vm.count("outfile")) BOOST_THROW_EXCEPTION(po::required_option("outfile"));
      if (!vm.count("outfreq")) BOOST_THROW_EXCEPTION(po::required_option("outfreq"));
      outfile = vm["outfile"].as<std::string>();
      outfreq = vm["outfreq"].as<int>();
    }

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
      run<kin_cloud_2d_blk_1m<setup::real_t, n_iters, ix>>(nx, nz, nt, outfile, outfreq);
    }
    else
    if (micro == "blk_2m")
    {
      struct ix { enum {rhod_th, rhod_rv, rhod_rc, rhod_rr, rhod_nc, rhod_nr}; };
      run<kin_cloud_2d_blk_2m<setup::real_t, n_iters, ix>>(nx, nz, nt, outfile, outfreq);
    }
    else 
    if (micro == "lgrngn")
    {
      struct ix { enum {rhod_th, rhod_rv}; };
      run<kin_cloud_2d_lgrngn<setup::real_t, n_iters, ix>>(nx, nz, nt, outfile, outfreq);
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
