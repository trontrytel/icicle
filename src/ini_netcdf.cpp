/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date January 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#include "cfg/cfg_netcdf.hpp"
#ifdef USE_NETCDF
#  include "ini_netcdf.hpp"
#  include "cmn/cmn_netcdf.hpp"
#  include "cmn/cmn_error.hpp"

#  include <vector>
using std::vector;

// pimpl struct definition
template <typename real_t>
struct ini_netcdf<real_t>::detail
{
  unique_ptr<NcFile> f;
  size_t xs, ys, zs; 
  string filename;
  const grd<real_t> *grid;
};

// ctor
template <typename real_t>
ini_netcdf<real_t>::ini_netcdf(const grd<real_t> &grid, const string filename)
  : pimpl(new detail())
{ 
  pimpl->filename = filename;
  pimpl->grid = &grid;
    try
    {
      pimpl->f.reset(new NcFile(filename, NcFile::read));
      // sanity checks if the grid matches current grid
      NcDim 
        d_x = pimpl->f->getDim("X"),
        d_y = pimpl->f->getDim("Y"),
        d_z = pimpl->f->getDim("Z"); // TODO: theta, r, phi
      if (d_x.isNull()) error_macro("X dimension not found in file " << filename)
      if (d_y.isNull()) error_macro("Y dimension not found in file " << filename) // TODO: should not be needed if ny == 1...
      if (d_z.isNull()) error_macro("Z dimension not found in file " << filename)
      pimpl->xs = d_x.getSize();
      pimpl->ys = d_y.getSize();
      pimpl->zs = d_z.getSize();
      if (
        (pimpl->xs != grid.nx() && pimpl->xs != 1) ||
        (pimpl->ys != grid.ny() && pimpl->ys != 1) ||
        (pimpl->zs != grid.nz() && pimpl->zs != 1)
      ) error_macro(
        "X,Y or Z dim extent (" << pimpl->xs << "," << pimpl->ys << "," << pimpl->zs << ") " <<
        "does not match grid extent (" << grid.nx() << "," << grid.ny() << "," << grid.nz() << ")")
      // TODO: check variable X,Y,Z (dx, Arakawa-C setting etc)
    }
    catch (NcException& e) error_macro(e.what())
}

template <typename real_t>
void ini_netcdf<real_t>::populate_scalar_field(
  const string &varname,
  const mtx::idx &ijk,
  mtx::arr<real_t> &data
) const
{
  try
  {
    NcVar v = pimpl->f->getVar(varname);
    if (v.isNull()) 
      error_macro("variable " << varname << " not found in file " << pimpl->filename)
    // TODO: check if it has (x,y,z) dimensions

    vector<size_t> startp(3), countp(3, 1);
    // reading a range or just one element if e.g. a profile is to be copied along a dimension
    countp[0] = (pimpl->xs == 1 ? 1 : ijk.ubound(mtx::i) - ijk.lbound(mtx::i) + 1);
    startp[0] = (pimpl->xs == 1 ? 0 : ijk.lbound(mtx::i));
    // due to presence of halos the data to be stored is not contiguous, 
    // hence looping over the two major ranks
    for (int k_int = ijk.lbound(mtx::k); k_int <= ijk.ubound(mtx::k); ++k_int) // loop over "outer" dimension
    { 
      startp[2] = (pimpl->zs == 1 ? 0 : k_int);
      for (int j_int = ijk.lbound(mtx::j); j_int <= ijk.ubound(mtx::j); ++j_int)
      { 
        startp[1] = (pimpl->ys == 1 ? 0 : j_int);
        if (pimpl->xs != 1)
        {
          assert(data(ijk.i, j_int, k_int).isStorageContiguous());
          v.getVar(startp, countp, data(ijk.i, j_int, k_int).dataFirst());
        }
        else
        {
          for (int i_int = ijk.lbound(mtx::i); i_int <= ijk.ubound(mtx::i); ++i_int)
            v.getVar(startp, countp, &data(i_int, j_int, k_int));
        }
      }   
    }  
  } 
  catch (NcException& e) error_macro(e.what())
}
#endif

// explicit instantiations
#include "cfg/cfg_types.hpp"
#if defined(USE_FLOAT)
template class ini_netcdf<float>;
#endif
#if defined(USE_DOUBLE)
template class ini_netcdf<double>;
#endif
#if defined(USE_LDOUBLE)
template class ini_netcdf<long double>;
#endif
