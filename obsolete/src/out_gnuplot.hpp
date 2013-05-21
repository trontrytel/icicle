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
#pragma once
#include "out.hpp"
#include "grd.hpp"

template <typename real_t>
class out_gnuplot : public out<real_t>
{
  private: unsigned long t_last;
  private: int i_last;
  private: bool steps;
  private: const grd<real_t> &grid;

  // ctor
  public: out_gnuplot(const grd<real_t> &grid, string usng);

  public: void record(
    const string &name,
    const mtx::arr<real_t> &psi,
    const mtx::idx &ijk, 
    const unsigned long t
  );
  
  private:
  template<class idx>
  void record_helper(const mtx::arr<real_t> &psi, const mtx::rng &i, const unsigned long t);
};
