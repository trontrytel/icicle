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
  private: const grd<real_t> *grid;
  public: ini_netcdf(const grd<real_t> &grid, const string filename)
    : grid(&grid)
  { 
    f.reset(new NcFile(filename, NcFile::read));
  }

  public: virtual void populate_scalar_field(
    const string &varname,
    const mtx::idx &ijk,
    mtx::arr<real_t> &data
  ) 
  {
cerr << "loading: " << varname << endl;
/*
    for (int i = ijk.lbound(mtx::i); i <= ijk.ubound(mtx::i); ++i)
      for (int j = ijk.lbound(mtx::j); j <= ijk.ubound(mtx::j); ++j)
        for (int k = ijk.lbound(mtx::k); k <= ijk.ubound(mtx::k); ++k)
          data(i,j,k) = psi(grid->x(i,j,k), grid->y(i,j,k), grid->z(i,j,k));
*/
  }
};
#endif
