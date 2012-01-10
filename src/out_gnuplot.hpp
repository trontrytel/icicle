/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
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
  private: unsigned long t_last; // TODO: si::seconds?
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
    arr<real_t> *psi,
    const rng &i, const rng &j, const rng &k, const unsigned long t
  ) 
  {
    // sanity check
    if (j.first() == j.last() && k.first() == k.last()) 
      record_helper<idx_ijk>(psi, i, t);
    else if (i.first() == i.last() && k.first() == k.last()) 
      record_helper<idx_jki>(psi, j, t);
    else if (i.first() == i.last() && j.first() == j.last()) 
      record_helper<idx_kij>(psi, k, t);
    else 
      error_macro("gnuplot output works for 1D simulations only") 
  }
  
  private:
  template<class idx>
  void record_helper(arr<real_t> *psi, const rng &i, const unsigned long t)
  {
    // end record if needed and output the data + some housekeeping
    if (t_last != -1 && t != t_last)
    {
      cout << endl << endl;
      i_last = -1;
    }
    assert(i_last + 1 == i.first()); // if data is output in order
    for (int ii=i.first(); ii<=i.last(); ++ii) 
    {
      if (steps)
        cout 
          << real_t((grid->x(ii,0,0) - real_t(.5) * grid->dx()) / si::metres) // TODO: this assumes x 
          << "\t"
          << *(*psi)(idx(rng(ii,ii), rng(0,0), rng(0,0))).dataFirst() 
          << endl
          << real_t((grid->x(ii,0,0) + real_t(.5) * grid->dx()) / si::metres) // TODO: this assumes x
          << "\t"
          << *(*psi)(idx(rng(ii,ii), rng(0,0), rng(0,0))).dataFirst() 
          << endl;
      else
        cout << *(*psi)(idx(rng(ii,ii), rng(0,0), rng(0,0))).dataFirst() << endl;
    }
  
    // housekeeping
    t_last = t;
    i_last = i.last();
  }
};

#endif
