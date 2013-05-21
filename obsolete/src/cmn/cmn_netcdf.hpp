/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date August 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @brief netcdf includes
 */

#pragma once

// the Lynton Appel's netCDF-4 C++ API (since netCDF 4.1.1)
#    ifdef USE_BOOST_MPI 
// fixes preprocessor macro redefinition conflict with MPI
// cf. http://www.unidata.ucar.edu/mailing_lists/archives/netcdfgroup/2009/msg00350.html
#      include "mpi/mpi.h"
#      define MPI_INCLUDED
#    endif    
#    include <netcdf>
using netCDF::NcFile;
using netCDF::NcVar;
using netCDF::NcDim;
using netCDF::ncFloat;
using netCDF::exceptions::NcException;
