/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 *  @section DESCRIPTION
 *    provides gnuplot-digestible ascii-based output facility suitable
 *    for quick-looking the results in 1D, e.g.:
 *    gnuplot> splot '< ./icicle --nx 100 --ny 1 --nz 1 ... ' u 0:-2:1 w lines palette
 */
#ifndef OUT_GNUPLOT_HPP
#  define OUT_GNUPLOT_HPP

#  include "out.hpp"

template <typename real_t>
class out_gnuplot : public out<real_t>
{
  private: unsigned long t_last; // TODO: si::seconds?
  private: int i_last;

  public: out_gnuplot() 
  { 
    t_last = -1; 
    i_last = -1; 
  }

  public: virtual void record(
    Array<real_t, 3> *psi,
    const Range &i, const Range &j, const Range &k, const unsigned long t
  ) 
  {
    // sanity check
    if (j.first() != j.last() || k.first() != k.last())
      error_macro("gnuplot output works for 1D simulations only") 

    // end record if needed and output the data + some housekeeping
    if (t_last != -1 && t != t_last)
    {
      cout << endl << endl;
      i_last = -1;
    }
    assert(i_last + 1 == i.first()); // if data is output in order
    Array<real_t, 1> a = (*psi)(i, 0, 0);
    for (int ii=a.lbound(0); ii<=a.ubound(0); ++ii) cout << a(ii) << endl;
  
    // housekeeping
    t_last = t;
    i_last = i.last();
  }
};

#endif
