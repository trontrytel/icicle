/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#include "out_gnuplot.hpp"
#include "cmn/cmn_error.hpp"

template <typename real_t> 
out_gnuplot<real_t>::out_gnuplot(const grd<real_t> &grid, string usng) 
  : t_last(-1), i_last(-1), grid(grid)
{ 
    if (usng == "1:-2:2") steps = true;
    else if (usng == "0:-2:1") steps = false;
    else error_macro("unknown using specification: " << usng << "(should be 1:-2:2 or 0:-2:1)")
}

template <typename real_t>
void out_gnuplot<real_t>::record(
    const string &name,
    const mtx::arr<real_t> &psi,
    const mtx::idx &ijk, 
    const unsigned long t
) 
{
    assert(name == "psi"); // TODO: maybe allow chosing which to output, or try outputting multiple fields?
    // sanity check
    if (ijk.lbound(mtx::j) == ijk.ubound(mtx::j) && ijk.lbound(mtx::k) == ijk.ubound(mtx::k)) 
      record_helper<mtx::idx_ijk>(psi, ijk.i, t);
    else if (ijk.lbound(mtx::i) == ijk.ubound(mtx::i) && ijk.lbound(mtx::k) == ijk.ubound(mtx::k)) 
      record_helper<mtx::idx_jki>(psi, ijk.j, t);
    else if (ijk.lbound(mtx::i) == ijk.ubound(mtx::i) && ijk.lbound(mtx::j) == ijk.ubound(mtx::j)) 
      record_helper<mtx::idx_kij>(psi, ijk.k, t);
    else 
      error_macro("gnuplot output works for 1D simulations only") 
}
  
template <typename real_t>
template<class idx>
void out_gnuplot<real_t>::record_helper(const mtx::arr<real_t> &psi, const mtx::rng &i, const unsigned long t)
{
    // end record if needed and output the data + some housekeeping
    if (t_last != -1 && t != t_last)
    {
      std::cout << std::endl << std::endl;
      i_last = -1;
    }
    assert(i_last + 1 == i.first()); // if data is output in order
    for (int ii=i.first(); ii<=i.last(); ++ii) 
    {
      if (steps)
        std::cout 
          << real_t((grid.x(ii,0,0) - real_t(.5) * grid.dx()) / si::metres) // TODO: this assumes x 
          << "\t"
          << *(psi)(idx(mtx::rng(ii,ii), mtx::rng(0,0), mtx::rng(0,0))).dataFirst() 
          << std::endl
          << real_t((grid.x(ii,0,0) + real_t(.5) * grid.dx()) / si::metres) // TODO: this assumes x
          << "\t"
          << *(psi)(idx(mtx::rng(ii,ii), mtx::rng(0,0), mtx::rng(0,0))).dataFirst() 
          << std::endl;
      else
        std::cout << *(psi)(idx(mtx::rng(ii,ii), mtx::rng(0,0), mtx::rng(0,0))).dataFirst() << std::endl;
    }
  
    // housekeeping
    t_last = t;
    i_last = i.last();
}

// explicit instantiations
#define ICICLE_INSTANTIATE_CLASS out_gnuplot
#include "cmn/cmn_instant.hpp"
