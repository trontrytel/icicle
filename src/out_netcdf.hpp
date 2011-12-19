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

// fixes preprocessor macro redefinition conflict with MPI
// cf. http://www.unidata.ucar.edu/mailing_lists/archives/netcdfgroup/2009/msg00350.html
//#    ifdef USE_BOOST_MPI 
//#      define MPI_INCLUDED
//#    endif              

// the Lynton Appel's netCDF-4 C++ API (since netCDF 4.1.1)
#    include <netcdf>
using namespace netCDF;
using namespace netCDF::exceptions;

// TODO: #include <boost/timer/timer.hpp> (requires Boost 1.48)
// TODO: ncFloat vs. ncDouble, ... ?
// TODO: add X_sclr i X_vctr variables! (e.g. for axis labelling)
// TODO: is the order of dimensions optimal?

template <typename real_t>
class out_netcdf : public out<real_t>
{
  private: auto_ptr<NcFile> f;
  private: NcVar vpsi;
  private: int freq; 

  public: out_netcdf(string file, grd<real_t> *grid, int nx, int ny, int nz, int freq, int ver) 
    : freq(freq)
  { 
    try
    {
      netCDF::NcFile::FileFormat fmt;
      switch (ver)
      {
        case 3: fmt = NcFile::classic; break;
        case 4: fmt = NcFile::nc4; break;
        default: error_macro("unsupported netCDF format version: " << ver)
      }
      f.reset(new NcFile(file, NcFile::newFile, fmt)); 
      NcDim 
        d_t = f->addDim("time"),
        d_xs = f->addDim("X", nx),
        d_ys = f->addDim("Y", ny),
        d_zs = f->addDim("Z", nz);
      {
        vector<NcDim> sdims(4);
        sdims[0] = d_t;
        sdims[1] = d_xs;
        sdims[2] = d_ys;
        sdims[3] = d_zs;
        vpsi = f->addVar("psi", ncFloat, sdims); 
      }
      {
        vector<NcDim> vdims(3);
        //vdims[0] = d_xv;
        //vdims[1] = d_yv;
        //vdims[2] = d_zv;
        //NcVar vu = f->addVar("u", ncFloat, vdims); 
        //NcVar vv = f->addVar("v", ncFloat, vdims); 
        //NcVar vw = f->addVar("w", ncFloat, vdims); 
      }
      // workaround for the lack of netCDF3 nc_enddef() in the netCDF C++4 API
      if (fmt == NcFile::classic) 
      {
        f.reset(); // closing the file
        f.reset(new NcFile(file, NcFile::write, fmt)); // reopening it
        vpsi = f->getVar("psi");
      }
    }
    catch (NcException& e) error_macro(e.what())
  }

  public: virtual void record(
    arr<real_t> *psi,
    const Range &i, const Range &j, const Range &k, const unsigned long t
  ) 
  {
    if (t % freq != 0) return;
    vector<size_t> startp(4), countp(4, 1);
    startp[0] = t / freq;
    countp[3] = k.last() - k.first() + 1;
    // due to presence of halos the data to be stored is not contiguous, 
    // hence looping over the two major ranks
    for (int i_int = i.first(); i_int <= i.last(); ++i_int) // loop over "outer" dimension
    {
      startp[1] = i_int;
      for (int j_int = j.first(); j_int <= j.last(); ++j_int)
      {
        assert((*psi)(i_int, j_int, k).isStorageContiguous());
        startp[2] = j_int;
        startp[3] = k.first();
        try 
        {
          vpsi.putVar(startp, countp, (*psi)(i_int, j_int, k).dataFirst());
        }
        catch (NcException& e) error_macro(e.what());
      }
    }
    //if (!f->sync()) warning_macro("failed to synchronise netCDF file")
  }
};
#  endif
#endif 
