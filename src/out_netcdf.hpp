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
#    include "inf.hpp"
#    include "eqs.hpp"
#    include "stp.hpp"
#    include "grd_arakawa-c-lorenz.hpp"

// fixes preprocessor macro redefinition conflict with MPI
// cf. http://www.unidata.ucar.edu/mailing_lists/archives/netcdfgroup/2009/msg00350.html
//#    ifdef USE_BOOST_MPI 
//#      define MPI_INCLUDED
//#    endif              

// the Lynton Appel's netCDF-4 C++ API (since netCDF 4.1.1)
#    include <netcdf>
using netCDF::NcFile;
using netCDF::NcVar;
using netCDF::NcDim;
using netCDF::ncFloat;
using netCDF::exceptions::NcException;

// TODO: ncFloat vs. ncDouble as a command-line option (but it float computations and double requeste -> error)
// TODO: add X_sclr i X_vctr variables! (e.g. for axis labelling)
// TODO: is the order of dimensions optimal?

template <typename real_t>
class out_netcdf : public out<real_t>
{
  private: auto_ptr<NcFile> f;
  private: vector<NcVar> vars;
  private: inf info;

  public: out_netcdf(
    const string &file, 
    stp<real_t> *setup, 
    int ver, 
    const string &cmdline
  ) 
    : info(cmdline)
  { 
    grd_arakawa_c_lorenz<real_t> *grid = dynamic_cast<grd_arakawa_c_lorenz<real_t>*>(setup->grid);
    if (grid == NULL) error_macro("netCDF output supports only the Arakawa-C grid")
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
        d_xs = f->addDim("X", grid->nx()),
        d_ys = f->addDim("Y", grid->ny()),
        d_zs = f->addDim("Z", grid->nz());
      {
        // dimensions
        vector<NcDim> sdims(4);
        sdims[0] = d_t;
        sdims[1] = d_xs;
        sdims[2] = d_ys;
        sdims[3] = d_zs; // TODO: skip dimensions of size 1?

        // grid variables
        NcVar v_xs = f->addVar("X", ncFloat, d_xs);
        v_xs.putAtt("unit", "metres");
        {
          real_t *tmp = new real_t[grid->nx()];
          for (int i = 0; i < grid->nx(); ++i) 
            tmp[i] = grid->i2x(i) / si::metres;
          v_xs.putVar(tmp);
          delete[] tmp;
        }
        NcVar v_ys = f->addVar("Y", ncFloat, d_ys);
        v_ys.putAtt("unit", "metres");
        { // TODO: merge with the above
          real_t *tmp = new real_t[grid->ny()];
          for (int j = 0; j < grid->ny(); ++j) 
            tmp[j] = grid->j2y(j) / si::metres;
          v_ys.putVar(tmp);
          delete[] tmp;
        }
        NcVar v_zs = f->addVar("Z", ncFloat, d_zs);
        v_zs.putAtt("unit", "metres");
        { // TODO: merge with the above
          real_t *tmp = new real_t[grid->nz()];
          for (int k = 0; k < grid->nz(); ++k) 
            tmp[k] = grid->k2z(k) / si::metres;
          v_zs.putVar(tmp);
          delete[] tmp;
        }

        // scalar fields
        for (int v = 0; v < setup->equations->n_vars(); ++v)
        {
          vars.push_back(f->addVar(setup->equations->var_name(v), ncFloat, sdims)); 
          vars.at(v).putAtt("unit", setup->equations->var_unit(v));
        }

        // Courant field TODO

        // timesteps
        NcVar 
          v_dtadv = f->addVar("dt_adv", ncFloat, vector<NcDim>()),
          v_dtout = f->addVar("dt_out", ncFloat, vector<NcDim>());
        v_dtadv.putAtt("unit", "seconds");
        v_dtout.putAtt("unit", "seconds");
        {
          real_t tmp = setup->dt / si::seconds;
          v_dtadv.putVar(&tmp);
        }
        { 
          real_t tmp = setup->dt * real_t(setup->nout) / si::seconds;
          v_dtout.putVar(&tmp);
        }
      }
    }
    catch (NcException& e) error_macro(e.what())
  }

  public: ~out_netcdf()
  {
    // from a destructor results in undefined behaviour 
    // so trying to handle the exceptions here
    try
    {
      // TODO: distmem concurency logic! (e.g. MPI)
      map<string,string> im = info.get_map();
      for (map<string,string>::iterator it = im.begin(); it != im.end(); ++it) 
        f->putAtt(it->first, it->second);
    }
    catch (NcException& e) warning_macro(e.what())
  }

  public: virtual void record(
    int e,
    const mtx::arr<real_t> &psi,
    const mtx::idx &ijk, 
    const unsigned long t // t is the number of the record!
  ) 
  {
    vector<size_t> startp(4), countp(4, 1);
    startp[0] = t;
    countp[3] = ijk.ubound(mtx::k) - ijk.lbound(mtx::k) + 1;
    // due to presence of halos the data to be stored is not contiguous, 
    // hence looping over the two major ranks
    for (int i_int = ijk.lbound(mtx::i); i_int <= ijk.ubound(mtx::i); ++i_int) // loop over "outer" dimension
    {
      startp[1] = i_int;
      for (int j_int = ijk.lbound(mtx::j); j_int <= ijk.ubound(mtx::j); ++j_int)
      {
        assert((psi)(i_int, j_int, ijk.k).isStorageContiguous());
        startp[2] = j_int;
        startp[3] = ijk.lbound(mtx::k);
        try 
        {
          vars.at(e).putVar(startp, countp, (psi)(i_int, j_int, ijk.k).dataFirst());
        }
        catch (NcException& e) error_macro(e.what());
      }
    }
    //if (!f->sync()) warning_macro("failed to synchronise netCDF file")
  }
};
#  endif
#endif 
