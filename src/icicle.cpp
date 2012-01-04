/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
/** @mainpage notitle
 *  @section README README file (incl. installation instructions)
 *  @verbinclude "../README"
 *  @section HACKING HACKING file (coding conventions)
 *  @verbinclude "../HACKING"
 */

#include "cfg.hpp"
#include "cmn.hpp"
#define ICICLE_OPT_DESCS
#include "opt_adv.hpp"
#include "opt_grd.hpp"
#include "opt_ini.hpp"
#include "opt_out.hpp"
#include "opt_slv.hpp"
#include "opt_vel.hpp"
#undef ICICLE_OPT_DESCS

extern void mdl_flt(const po::variables_map&, const string&);
extern void mdl_dbl(const po::variables_map&, const string&);
extern void mdl_ldb(const po::variables_map&, const string&);

int main(int ac, char* av[])
{
  cerr << "-- init: icicle starting (built on " << __DATE__ << ")" << endl;
  try
  {
    // options list
    po::options_description desc("options");
    desc.add_options()
      ("help", "print this message")
      ("bits", po::value<int>()->default_value(32), "floating point bits: sizeof(float), sizeof(double), sizeof(long double)");
// <TODO> move somewhere else...
    desc.add_options()
      ("nsd", po::value<int>(), "number of subdomains (nx/nsd must be int)")
      ("nx", po::value<int>()->default_value(1), "number of grid points (X)")
      ("ny", po::value<int>()->default_value(1), "number of grid points (Y)")
      ("nz", po::value<int>()->default_value(1), "number of grid points (Z)")
      ("nt", po::value<unsigned long>(), "number of timesteps")
      ("dt", po::value<string>(), "timestep length [s]");
// </TODO> 
    opt_adv_desc(desc);
    opt_slv_desc(desc);
    opt_out_desc(desc);
    opt_grd_desc(desc);
    opt_vel_desc(desc);
    opt_ini_desc(desc);
    po::variables_map vm;
    po::store(po::parse_command_line(ac, av, desc), vm);
    po::notify(vm);

    cerr << "-- init: parsing command-line options" << endl;

    // --help or no argument case
    if (vm.count("help") || ac == 1)
    {
      cerr << desc << endl;
      exit(EXIT_FAILURE);
    }

    // --slv list
    if (vm.count("slv") && vm["slv"].as<string>() == "list")
    {
      cout << "serial fork";
#  ifdef _OPENMP
      cout << " openmp fork+openmp";
#  endif
#  ifdef USE_BOOST_THREAD
      cout << " threads fork+threads";
#  endif
#  ifdef USE_BOOST_MPI
      cout << " mpi";
#    ifdef USE_BOOST_THREAD
      cout << " mpi+threads";
#    endif
#    ifdef _OPENMP
      cout << " mpi+openmp";
#    endif
#  endif
      cout << endl;
      exit(EXIT_FAILURE);
    }

    // string containing all passed options (e.g. for archiving in a netCDF file)
    ostringstream options;
    options << string(av[0]);
    for (int i = 1; i < ac; ++i) options << string(" ") << string(av[i]);

    // --bits (floating point precision choice)
    int bits = vm["bits"].as<int>();
    if (sizeof(float) * 8 == bits) mdl_flt(vm, options.str());
    else if (sizeof(double) * 8 == bits) mdl_dbl(vm, options.str());
    else if (sizeof(long double) * 8 == bits) mdl_ldb(vm, options.str());
    else error_macro("unsupported number of bits (" << bits << ")")
  }
  catch (exception &e)
  {
    cerr << "-- exception cought: " << e.what() << endl;
    cerr << "-- exit: KO" << endl;
    exit(EXIT_FAILURE);
  }
  cerr << "-- exit: OK" << endl;
  exit(EXIT_SUCCESS);
}
