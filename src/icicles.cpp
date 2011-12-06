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
  cerr << "-- init: icicles starting (built on " << __DATE__;
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
      ("bits", po::value<int>()->default_value(32), "floating point bits: sizeof(float), sizeof(double), sizeof(long double)")
      ("adv", po::value<string>(), "advection scheme: leapfrog, upstream, mpdata")
      ("adv.mpdata.iord", po::value<int>(), "mpdata iord option: 1, 2, ...")
      ("slv", po::value<string>(), "solver: serial, openmp, threads")
      ("out", po::value<string>(), "output: debug, gnuplot, netcdf")
      ("out.netcdf.file", po::value<string>(), "output filename (e.g. for netcdf)")
      ("out.netcdf.freq", po::value<int>(), "inverse of output frequency (10 for every 10-th record)")
      ("out.netcdf.ver", po::value<int>()->default_value(4), "netCDF format version (3 or 4)")
      ("nsd", po::value<int>(), "number of subdomains (nx/nsd must be int)")
      ("nx", po::value<int>(), "number of grid points (X)")
      ("ny", po::value<int>()->default_value(1), "number of grid points (Y)")
      ("nz", po::value<int>()->default_value(1), "number of grid points (Z)")
      ("nt", po::value<unsigned long>(), "number of timesteps")
      ("dt", po::value<string>(), "timestep length [s]")
      ("grd", po::value<string>()->default_value("arakawa-c-lorenz"), "grid: arakawa-c-lorenz")
      ("grd.dx", po::value<string>(), "gridbox length (X) [m]")
      ("grd.dy", po::value<string>()->default_value("1"), "gridbox length (Y) [m]")
      ("grd.dz", po::value<string>()->default_value("1"), "gridbox length (Z) [m]")
      ("vel", po::value<string>(), "velocity field: uniform, rasinski")
      ("vel.uniform.u", po::value<string>(), "velocity (X) [m/s]")
      ("vel.uniform.v", po::value<string>()->default_value("0"), "velocity (Y) [m/s]")
      ("vel.uniform.w", po::value<string>()->default_value("0"), "velocity (Z) [m/s]")
      ("vel.rasinski.Z_clb", po::value<string>(), "cloud base height [m]")
      ("vel.rasinski.Z_top", po::value<string>(), "cloud top height [m]")
      ("vel.rasinski.A", po::value<string>(), "amplitude [m2/s]")
      ("vel.test.omega", po::value<string>(), "frequency [1/s]")
      ("ini", po::value<string>(), "initical condition: origin-one, cone")
      ("ini.cone.h", po::value<string>(), "h [m]")
      ("ini.cone.x0", po::value<string>(), "x0 [m]")
      ("ini.cone.z0", po::value<string>(), "z0 [m]")
      ("ini.cone.r", po::value<string>(), "r [m]")
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

    // --bits (floating point precision choice)
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
