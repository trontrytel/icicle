/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
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
  private: auto_ptr<NcFile> f;
  private: NcVar *vpsi;

  public: out_netcdf(string file, grd<real_t> *grid, int nx, int ny, int nz) 
  { 
    f.reset(new NcFile(file.c_str(), NcFile::New)); // TODO: other parameters (perhaps via variables_map?)
    if (!f->is_valid()) error_macro("failed to open netcdf file for writing: " << file)
    NcDim 
      *t = f->add_dim("T"),
      *xs = f->add_dim("X_sclr", nx),
      *ys = f->add_dim("Y_sclr", ny),
      *zs = f->add_dim("Z_sclr", nz),
      *xv = f->add_dim("X_vctr", grid->rng_vctr(0, nx-1).length()), // TODO:  
      *yv = f->add_dim("Y_vctr", grid->rng_vctr(0, ny-1).length()), // TODO: that's a kludge
      *zv = f->add_dim("Z_vctr", grid->rng_vctr(0, nz-1).length()); // TODO:
    // TODO: is the order of dimensions optimal?
    vpsi = f->add_var("psi", ncFloat, t, xs, ys, zs); // TODO: ncFloat vs. ncDouble, ...
    NcVar *vu = f->add_var("u", ncFloat, xv, yv, zv); // TODO: ncFloat vs. ncDouble, ...
    NcVar *vv = f->add_var("v", ncFloat, xv, yv, zv); // TODO: ncFloat vs. ncDouble, ...
    NcVar *vw = f->add_var("w", ncFloat, xv, yv, zv); // TODO: ncFloat vs. ncDouble, ...
// TODO: add X_sclr i X_vctr variables! (e.g. for axis labelling)

    // a sanity check to verify if Boost.units was optimised correctly and if pointer
    // arithmetics may be applied to &(blitz::Array<boost::units::quantity>(...).value())
    {
      // TODO: is it really needed???
      Array<quantity<unit, real_t>, 1> a(3);
      a = 11,22,33;
      if (*(&a(0).value() + 2) != 33) error_macro("The compiler did not optimise Blitz+Boost.Units enough :(");
    }
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
