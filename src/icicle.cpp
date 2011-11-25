/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */

#include "config.hpp"
#include "common.hpp"
#include "model.hpp"

int main(int ac, char* av[])
{
  cerr << "-- init: icicle starting (built on " << __DATE__;
#ifdef __FAST_MATH__
  cerr << " with FAST_MATH enabled!";
#endif
  cerr << ")" << endl;
  try
  {
    // options list
    po::options_description desc("options");
    desc.add_options()
      ("help", "print this message")
      ("bits", po::value<int>(), "floating point bits: sizeof(float), sizeof(double), sizeof(long double)")
      ("adv", po::value<string>(), "advection scheme: leapfrog, mpdata")
      ("adv.mpdata.iord", po::value<int>(), "mpdata iord option: 1, 2, ...")
      ("dom", po::value<string>(), "domain: serial, openmp, threads")
      ("out", po::value<string>(), "output: gnuplot, netcdf")
      ("out.netcdf.file", po::value<string>(), "output filename (e.g. for netcdf)")
      ("nsd", po::value<int>(), "number of subdomains (nx/nsd must be int)")
      ("nx", po::value<int>(), "number of grid points (X)")
      ("ny", po::value<int>(), "number of grid points (Y)")
      ("nz", po::value<int>(), "number of grid points (Z)")
      ("nt", po::value<unsigned long>(), "number of timesteps")
      ("dt", po::value<string>(), "timestep length [s]")
      ("grd", po::value<string>(), "grid: arakawa-c")
      ("grd.dx", po::value<string>(), "gridbox length (X) [m]")
      ("grd.dy", po::value<string>(), "gridbox length (Y) [m]")
      ("grd.dz", po::value<string>(), "gridbox length (Z) [m]")
      ("vel", po::value<string>(), "velocity field: uniform, rasinski")
      ("vel.uniform.u", po::value<string>(), "velocity (X) [m/s]")
      ("vel.uniform.v", po::value<string>(), "velocity (Y) [m/s]")
      ("vel.uniform.w", po::value<string>(), "velocity (Z) [m/s]")
      ("vel.rasinski.Z_clb", po::value<string>(), "cloud base height [m]")
      ("vel.rasinski.Z_top", po::value<string>(), "cloud top height [m]")
      ("vel.rasinski.A", po::value<string>(), "amplitude [m2/s]")
    ;
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

    // --dom list
    if (vm.count("dom") && vm["dom"].as<string>() == "list")
    {
      cout << "serial";// FIXME: fork";
#  ifdef _OPENMP
      cout << " openmp";
#  endif
#  ifdef USE_BOOST_THREADS
      cout << " threads";
#  endif
      cout << endl;
      exit(EXIT_FAILURE);
    }

    // --bits (floating point precision choice)
    if (!vm.count("bits")) error_macro("floating point precision not specified (--bits option)")
    int bits = vm["bits"].as<int>();
    if (sizeof(float) * 8 == bits) model<float>(vm);
//    else if (sizeof(double) * 8 == bits) model<double>(vm);
//    else if (sizeof(long double) * 8 == bits) model<long double>(vm);
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
