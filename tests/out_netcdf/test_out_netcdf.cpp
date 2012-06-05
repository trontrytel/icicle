/** @file
 *  @example out_netcdf/test_out_netcdf.cpp
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date December 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @section DESCRIPTION
 *    
 *
 *    Here's an example header of the netCDF file (output of ncdump -c):
 *    \include "out_netcdf/lst1.txt"
 */

#include <cstdlib>
#include <string>
extern "C" {
#include <unistd.h>
}

int main() 
{ 
  unlink("test_netcdf.nc");
  int status;
  std::string cmd =
    "../../icicle"
    " --grd.dx 1 --grd.dy 1 --grd.dz 1"
    " --vel uniform --vel.uniform.u -1"
    " --t_max 10 --dt_out 1 --grd.nx 20"
    " --adv mpdata"
    " --slv serial"
    " --out netcdf --out.netcdf.file test_netcdf.nc"
    " --ini boxcar --ini.boxcar.bx 1";
  if (0 != system(cmd.c_str())) 
    exit(EXIT_FAILURE);
  if (0 != system("test \"`ncdump -k test_netcdf.nc`\" = \"netCDF-4\"")) 
    exit(EXIT_FAILURE);
  if (0 != system("ncdump -c test_netcdf.nc > lst1.txt")) 
    exit(EXIT_FAILURE);
  unlink("test_netcdf.nc");

  if (0 != system((cmd + " --out.netcdf.ver 3").c_str())) 
    exit(EXIT_FAILURE);
  if (0 != system("test \"`ncdump -k test_netcdf.nc`\" = \"classic\"")) 
    exit(EXIT_FAILURE);
  unlink("test_netcdf.nc");
}
