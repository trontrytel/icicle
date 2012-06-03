/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @brief definition of and helper macros for \ref adv - the base class for advection operators
 */
#ifndef ADV_HPP
#  define ADV_HPP

#  include "cmn.hpp"
#  include "mtx.hpp"

/// @brief a base class for advection operators
template <typename real_t>
class adv 
{
  public: virtual const int stencil_extent() = 0;
  public: virtual const int time_levels() = 0;
  public: virtual const int num_steps() { return 1; }
  public: virtual const int num_sclr_caches() { return 0; }
  public: virtual const int num_vctr_caches() { return 0; }

  public: virtual const real_t courant_max() = 0; 
  public: virtual const real_t courant_min() = 0;

  // functor factory (TODO: rename to op3D)
  public: class op3D 
  {
    private: bool dox, doy, doz;
    public: op3D(const mtx::idx &ijk) :
      dox(ijk.i.first() != ijk.i.last()),
      doy(ijk.j.first() != ijk.j.last()),
      doz(ijk.k.first() != ijk.k.last())
    {}
    protected: bool do_x() { return dox; }
    protected: bool do_y() { return doy; }
    protected: bool do_z() { return doz; }
    public: virtual void operator()(
      mtx::arr<real_t> *psi[], 
      const int n, 
      const int s,
      const mtx::arr<real_t> * const Cx, 
      const mtx::arr<real_t> * const Cy, 
      const mtx::arr<real_t> * const Cz
    ) = 0;
  };

  public: virtual op3D *factory(
    const mtx::idx &ijk,
    mtx::arr<real_t> *tmp_s[], 
    mtx::arr<real_t> *tmp_v[],
    bool positive_definite
  ) = 0;

    // TODO: document, include as an option / autotest in CMake
#  ifdef MPDATA_NEGPOSPART_ABS
    protected: mtx_expr_1arg_macro(pospart, x,
      real_t(.5) * (x + abs(x))
    )   
    protected: mtx_expr_1arg_macro(negpart, x,
      real_t(.5) * (x - abs(x))
    )   
#  else
    protected: mtx_expr_1arg_macro(pospart, x,
      max(real_t(0), x)
    )   
    protected: mtx_expr_1arg_macro(negpart, x,
      min(real_t(0), x)
    )   
#  endif

};
#endif
