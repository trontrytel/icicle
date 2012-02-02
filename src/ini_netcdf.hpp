/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date January 2012
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef INI_NETCDF_HPP
#  define INI_NETCDF_HPP

#  include "ini.hpp"

template <typename real_t>
class ini_netcdf : public ini<real_t>
{
  private: auto_ptr<NcFile> f;
  private: string filename;
  private: const grd<real_t> *grid;
  public: ini_netcdf(const grd<real_t> &grid, const string filename)
    : filename(filename), grid(&grid)
  { 
    try
    {
      f.reset(new NcFile(filename, NcFile::read));
      // sanity checks if the grid matches current grid
      NcDim 
        d_x = f->getDim("X"),
        d_y = f->getDim("Y"),
        d_z = f->getDim("Z"); // TODO: theta, r, phi
      if (d_x.isNull()) error_macro("X dimension not found in file " << filename)
      if (d_y.isNull()) error_macro("Y dimension not found in file " << filename) // TODO: should not be needed if ny == 1...
      if (d_z.isNull()) error_macro("Z dimension not found in file " << filename)
      if (
        d_x.getSize() != grid.nx() ||
        d_y.getSize() != grid.ny() ||
        d_z.getSize() != grid.nz()
      ) error_macro(
        "X,Y or Z dim extent (" << d_x.getSize() << "," << d_y.getSize() << "," << d_z.getSize() << ") " <<
        "does not match grid extent (" << grid.nx() << "," << grid.ny() << "," << grid.nz() << ")")
      // TODO: check variable X,Y,Z (dx, Arakawa-C setting etc)
    }
    catch (NcException& e) error_macro(e.what())
  }

  public: virtual void populate_scalar_field(
    const string &varname,
    const mtx::idx &ijk,
    mtx::arr<real_t> &data
  ) 
  {
    try
    {
      NcVar v = f->getVar(varname);
      if (v.isNull()) 
        error_macro("variable " << varname << " not found in file " << filename)
      // TODO: check if it has (x,y,z) dimensions
      vector<size_t> startp(3), countp(3, 1);
      countp[2] = ijk.ubound(mtx::k) - ijk.lbound(mtx::k) + 1;
      // due to presence of halos the data to be stored is not contiguous, 
      // hence looping over the two major ranks
      for (int i_int = ijk.lbound(mtx::i); i_int <= ijk.ubound(mtx::i); ++i_int) // loop over "outer" dimension
      { 
        startp[0] = i_int;
        for (int j_int = ijk.lbound(mtx::j); j_int <= ijk.ubound(mtx::j); ++j_int)
        { 
          assert(data(i_int, j_int, ijk.k).isStorageContiguous());
          startp[1] = j_int;
          startp[2] = ijk.lbound(mtx::k);
          v.getVar(startp, countp, data(i_int, j_int, ijk.k).dataFirst());
        }   
      }  
    } 
    catch (NcException& e) error_macro(e.what())
  }
};
#endif
