/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 @  @brief contains the @ref main() function in which the floating point precision choice takes place
 */

#include "mdl.hpp"
#include "cfg/cfg_types.hpp"
#include "cfg/cfg_boost_mpi.hpp"
#include "cfg/cfg_boost_thread.hpp"

int main(int ac, char* av[])
{
  notice_macro("icicle starting (built on " << __DATE__ << ")")
#if defined(__GNUC__) && !defined(__FAST_MATH__)
  warning_macro("GCC was used without the -ffast-math flag!")
#endif
  try
  {
    // options list
    po::options_description desc("options");
    desc.add_options()
      ("help", "print this message")
      ("bits", po::value<int>()->default_value(
#if defined(USE_DOUBLE)
        8 * sizeof(double)
#elif defined(USE_FLOAT)
        8 * sizeof(float)
#elif defined(USE_LDOUBLE)
        8 * sizeof(long double)
#elif defined(USE_FLOAT128)
        8 * sizeof(__float128)
#else
#  error
#endif
      ), "floating point bits: sizeof(float), sizeof(double), sizeof(long double)");
    opt_stp_desc(desc);
    opt_adv_desc(desc);
    opt_slv_desc(desc);
    opt_out_desc(desc);
    opt_grd_desc(desc);
    opt_vel_desc(desc);
    opt_ini_desc(desc);
    opt_eqs_desc(desc);
    po::variables_map vm;
    po::store(po::parse_command_line(ac, av, desc), vm);
    po::notify(vm);

    // --help or no argument case
    if (vm.count("help") || ac == 1)
    {
      std::cerr << desc << std::endl;
      exit(EXIT_SUCCESS); // this is what GNU coding standards suggest
    }

    // string containing all passed options (e.g. for archiving in a netCDF file)
    ostringstream options;
    options << string(av[0]);
    for (int i = 1; i < ac; ++i) options << string(" ") << string(av[i]);

    // --bits (floating point precision choice)
    int bits = vm["bits"].as<int>();
#ifdef USE_FLOAT
    if (sizeof(float) * 8 == bits) 
      mdl<float>(vm, options.str());
    else 
#endif
#ifdef USE_DOUBLE
    if (sizeof(double) * 8 == bits) 
      mdl<double>(vm, options.str());
    else 
#endif
#ifdef USE_LDOUBLE
    if (sizeof(long double) * 8 == bits) 
      mdl<long double>(vm, options.str());
    else 
#endif
#ifdef USE_FLOAT128
    if (sizeof(__float128) * 8 == bits) 
      mdl<__float_128>(vm, options.str());
    else 
#endif
    error_macro("unsupported number of bits (" << bits << ")")
  }
  catch (std::exception &e)
  {
    std::cerr << "-- exception cought: " << e.what() << std::endl;
    std::cerr << "-- exit: KO" << std::endl;
    exit(EXIT_FAILURE);
  }
  std::cerr << "-- exit: OK" << std::endl;
  exit(EXIT_SUCCESS);
}
