/** @file
 *  @author Sylwester Arabas (University of Warsaw) <slayoo@igf.fuw.edu.pl>
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef OUT_NETCDF_HPP
#  define OUT_NETCDF_HPP
#  ifdef USE_NETCDF

#    include "out.hpp"
#    include <netcdfcpp.h>

// TODO: error handling

template <class unit, typename real_t>
class out_netcdf : public out<unit, real_t>
{
  private: NcFile *f;
  private: NcVar *vpsi;

  public: out_netcdf(string file, int nx, int ny, int nz) 
  { 
    f = new NcFile(file.c_str(), NcFile::New); // TODO: other parameters (perhaps via variables_map?)
    if (!f->is_valid()) error_macro("failed to open netcdf file for writing: " << file)
    NcDim 
      *t = f->add_dim("T"),
      *x = f->add_dim("X", nx),
      *y = f->add_dim("Y", ny),
      *z = f->add_dim("Z", nz);
    vpsi = f->add_var("psi", ncFloat, t, x, y, z); // TODO: ncFloat vs. ncDouble, ...

    // a sanity check to verify if Boost.units was optimised correctly and if pointer
    // arithmetics may be applied to &(blitz::Array<boost::units::quantity>(...).value())
    {
      // TODO: is it really needed???
      Array<quantity<unit, real_t>, 1> a(3);
      a = 11,22,33;
      if (*(&a(0).value() + 2) != 33) error_macro("The compiler did not optimise Blitz+Boost.Units enough :(");
    }
  }

  public: ~out_netcdf()
  {
    delete f;
  }

  public: virtual void record(
    Array<quantity<unit, real_t>, 3> **psi, const int n, 
    const Range &i, const Range &j, const Range &k, const unsigned long t
  ) 
  {
    // due to presence of halos the data to be stored is not contiguous, 
    // hence looping over the two major ranks
    for (int i_int = i.first(); i_int <= i.last(); ++i_int) // loop over "outer" dimension
    {
      for (int j_int = j.first(); j_int <= j.last(); ++j_int)
      {
        assert((*psi[n])(i_int, j_int, k).isStorageContiguous());
        if (!vpsi->set_cur(t, i_int, j_int, k.first()))
          error_macro("failed to set position in the netCDF file")
        if (!vpsi->put(
          &(*psi[n])(i_int, j_int, k).dataFirst()->value(), 1, 1, 1, (k.last() - k.first() + 1)
        )) error_macro("failed to write to netCDF file");
      }
    }
    if (!f->sync()) warning_macro("failed to synchronise netCDF file")
  }
};
#  endif
#endif 
