/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @section DESCRIPTION
 *    provides gnuplot-digestible ascii-based output facility suitable
 *    for quick-looking the results in 1D, e.g.:
 *    gnuplot> splot '< ./icicle --out.gnuplot' u 0:-2:1 w lines palette
 *    gnuplot> splot '< ./icicle --out.gnuplot --out.gnuplot.using 0:-2:1' u 0:-2:1 w lines palette
 *    gnuplot> splot '< ./icicle --out.gnuplot --out.gnuplot.using 1:-2:2' u 1:-2:2 w lines palette
 */
#ifndef OUT_GNUPLOT_HPP
#  define OUT_GNUPLOT_HPP

#  include "out.hpp"
#  include "grd.hpp"

template <typename real_t>
class out_gnuplot : public out<real_t>
{
  private: unsigned long t_last;
  private: int i_last;
  private: bool steps;
  private: grd<real_t> *grid;

  public: out_gnuplot(grd<real_t> *grid, string usng) 
    : t_last(-1), i_last(-1), grid(grid)
  { 
    if (usng == "1:-2:2") steps = true;
    else if (usng == "0:-2:1") steps = false;
    else error_macro("unknown using specification: " << usng << "(should be 1:-2:2 or 0:-2:1)")
  }

  public: virtual void record(
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
  
  private:
  template<class idx>
  void record_helper(const mtx::arr<real_t> &psi, const mtx::rng &i, const unsigned long t)
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
          << real_t((grid->x(ii,0,0) - real_t(.5) * grid->dx()) / si::metres) // TODO: this assumes x 
          << "\t"
          << *(psi)(idx(mtx::rng(ii,ii), mtx::rng(0,0), mtx::rng(0,0))).dataFirst() 
          << std::endl
          << real_t((grid->x(ii,0,0) + real_t(.5) * grid->dx()) / si::metres) // TODO: this assumes x
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
};

#endif
