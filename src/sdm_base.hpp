/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date March-June 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef SDM_BASE_HPP
#  define SDM_BASE_HPP

#  if defined(USE_THRUST) 
#    include <thrust/device_vector.h>
#    include <thrust/random.h> // TODO: random is not used in this file...
#    include <thrust/sequence.h>
#    include <thrust/sort.h>
#    include <thrust/iterator/constant_iterator.h>
#  endif

namespace sdm 
{

  typedef double thrust_real_t; // TODO: option / check if the device supports it
  typedef thrust::device_vector<int>::size_type thrust_size_t;

  // nested structure: super-droplet environment (velocity, temperature and moisture field)
  template <typename real_t>
  struct envi_t
  {
    // velocity field (copied from an Arakawa-C grid)
    thrust::device_vector<thrust_real_t> vx, vy;
    int vx_nx, vx_ny, vy_nx, vy_ny, n_cell; 

    // temperature and water vapour density fields (copied from an Arakawa-C grid)
    thrust::device_vector<thrust_real_t> rhod, rhod_th, rhod_rv;

    // ctor
    envi_t(int nx, int ny) 
    {
      vx_nx = nx + 1;
      vx_ny = ny;
      vy_nx = nx;
      vy_ny = ny + 1;
      vx.resize(vx_nx * vx_ny);
      vy.resize(vy_nx * vy_ny);

      n_cell = nx * ny;
      // TODO + 2 x halo if interpolating?
      rhod.resize(n_cell); 
      rhod_th.resize(n_cell); 
      rhod_rv.resize(n_cell);
    }
  };

  // nested structure: super-droplet state info
  template <typename real_t>
  struct stat_t
  {
    // number of particles (i.e. super-droplets)
    thrust_size_t n_part;

    // SD parameters that are variable from the ODE point of view (x,y,xi)
    thrust::device_vector<thrust_real_t> xy, xi; 

    // ... and since x and y are hidden in one SD.xy, we declare helper iterators
    thrust::device_vector<thrust_real_t>::iterator x_begin, x_end, y_begin, y_end;

    // SD parameters that are constant from the ODE point of view (n,rd3,i,j)
    thrust::device_vector<thrust_real_t> rd3, kpa; // rd3 -> dry radius to the third power 
    thrust::device_vector<thrust_size_t> id, n; // n -> number of real droplets in a super-droplet
    thrust::device_vector<int> i, j, ij; // location of super-droplet within the grid

    // sdm::stat ctor
    stat_t(int nx, int ny, real_t sd_conc_mean)
    {
      n_part = thrust_size_t(real_t(nx * ny) * sd_conc_mean);
      xy.resize(2 * n_part);
      xi.resize(n_part);
      rd3.resize(n_part);
      kpa.resize(n_part);
      i.resize(n_part);
      j.resize(n_part);
      ij.resize(n_part);
      n.resize(n_part);
      id.resize(n_part);

      x_begin = xy.begin(); 
      x_end   = x_begin + n_part;
      y_begin = x_end;
      y_end   = y_begin + n_part;
    }
  };

}
#endif
