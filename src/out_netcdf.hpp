/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef OUT_NETCDF_HPP
#  define OUT_NETCDF_HPP
#  ifdef USE_NETCDF

#    include "out.hpp"
#    include "inf.hpp"
#    include "eqs.hpp"
#    include "stp.hpp"
#    include "grd_arakawa-c-lorenz.hpp"

// TODO: ncFloat vs. ncDouble as a command-line option (but it float computations and double requeste -> error)
// TODO: add X_sclr i X_vctr variables! (e.g. for axis labelling)
// TODO: is the order of dimensions optimal?

template <typename real_t>
class out_netcdf : public out<real_t>
{
  private: unique_ptr<NcFile> f;
  private: map<string, NcVar> vars;
  private: inf info;

  public: out_netcdf(
    const string &file, 
    const stp<real_t> &setup, 
    int ver, 
    const string &cmdline
  ) 
    : info(cmdline)
  { 
    // TODO: that's about cartesian/spherical/etc, not about Arakawa-C
    try
    {
      netCDF::NcFile::FileFormat fmt;
      switch (ver)
      {
        case 3: fmt = NcFile::classic; break; // TODO: check in Cmake if that's supported!
        case 4: fmt = NcFile::nc4; break;
        default: error_macro("unsupported netCDF format version: " << ver)
      }
      f.reset(new NcFile(file, NcFile::newFile, fmt)); 

      NcDim 
        d_t = f->addDim("time"),
        d_xs = f->addDim("X", setup.grid.nx()),
        d_ys = f->addDim("Y", setup.grid.ny()),
        d_zs = f->addDim("Z", setup.grid.nz());
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
          real_t *tmp = new real_t[setup.grid.nx()];
          for (int i = 0; i < setup.grid.nx(); ++i) 
            tmp[i] = setup.grid.i2x(i) / si::metres;
          v_xs.putVar(tmp);
          delete[] tmp;
        }
        NcVar v_ys = f->addVar("Y", ncFloat, d_ys);
        v_ys.putAtt("unit", "metres");
        { // TODO: merge with the above
          real_t *tmp = new real_t[setup.grid.ny()];
          for (int j = 0; j < setup.grid.ny(); ++j) 
            tmp[j] = setup.grid.j2y(j) / si::metres;
          v_ys.putVar(tmp);
          delete[] tmp;
        }
        NcVar v_zs = f->addVar("Z", ncFloat, d_zs);
        v_zs.putAtt("unit", "metres");
        { // TODO: merge with the above
          real_t *tmp = new real_t[setup.grid.nz()];
          for (int k = 0; k < setup.grid.nz(); ++k) 
            tmp[k] = setup.grid.k2z(k) / si::metres;
          v_zs.putVar(tmp);
          delete[] tmp;
        }

        // scalar fields
        for (int v = 0; v < setup.eqsys.n_vars(); ++v)
        {
          string name = setup.eqsys.var_name(v);
          vars[name] = f->addVar(name, ncFloat, sdims); 
          vars[name].putAtt("unit", setup.eqsys.var_unit(v));
          vars[name].putAtt("description", setup.eqsys.var_desc(v));
        }

        // auxiliary fields
        for (const string &name : setup.eqsys.aux_names())
        {
          if (setup.eqsys.aux_tobeoutput(name)) 
          {
            // TODO: assert unique!
            vars[name] = f->addVar(name, ncFloat, sdims); 
            vars[name].putAtt("unit", setup.eqsys.aux_unit(name));
            vars[name].putAtt("description", setup.eqsys.aux_desc(name));
          }
        }

        // Courant field TODO

        // timesteps
        NcVar 
          v_dtadv = f->addVar("dt_adv", ncFloat, /* TODO: remove after upstream release */ vector<NcDim>()),
          v_dtout = f->addVar("dt_out", ncFloat, /* TODO: remove after upstream release */ vector<NcDim>());
        v_dtadv.putAtt("unit", "seconds");
        v_dtout.putAtt("unit", "seconds");
        {
          real_t tmp = setup.dt / si::seconds;
          v_dtadv.putVar(&tmp);
        }
        { 
          real_t tmp = setup.dt * real_t(setup.nout) / si::seconds;
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
    const string &name, 
    const mtx::arr<real_t> &psi,
    const mtx::idx &ijk, 
    const unsigned long t // t is the number of the record!
  ) 
  {
    try 
    {
      vector<size_t> startp(4), countp(4, 1);
      startp[0] = t;
      countp[1] = ijk.ubound(mtx::i) - ijk.lbound(mtx::i) + 1;
      startp[1] = ijk.lbound(mtx::i);
      // due to presence of halos the data to be stored is not contiguous, 
      // hence looping over the two major ranks
      for (int k_int = ijk.lbound(mtx::k); k_int <= ijk.ubound(mtx::k); ++k_int) // loop over "outer" dimension
      {
        startp[3] = k_int;
        for (int j_int = ijk.lbound(mtx::j); j_int <= ijk.ubound(mtx::j); ++j_int)
        {
          startp[2] = j_int;
          assert((psi)(ijk.i, j_int, k_int).isStorageContiguous());
          vars[name].putVar(startp, countp, (psi)(ijk.i, j_int, k_int).dataFirst());
        }
      }
      //if (!f->sync()) warning_macro("failed to synchronise netCDF file")
    }
    catch (NcException& e) error_macro(e.what());
  }
};
#  endif
#endif 
