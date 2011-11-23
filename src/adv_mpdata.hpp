/** @file
 *  @author Sylwester Arabas (University of Warsaw) <slayoo@igf.fuw.edu.pl>
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef ADV_MPDATA_HPP
#  define ADV_MPDATA_HPP

#  include "adv.hpp"
#  include "grd_2d-xz_arakawa-c.hpp"

template <class unit, typename real_t> 
class adv_mpdata : public adv<unit, real_t> 
{
  public: const int stencil_extent() { return 3; }
  public: const int time_levels() { return 2; }
  public: const int num_steps() { return iord; }

  private: int iord;
  private: grd_2d_xz_arakawa_c<real_t> *grid;

  public: adv_mpdata(grd<real_t> *g, int iord) 
    : iord(iord)
  {
    if (iord <= 0) error_macro("iord (the number of iterations) must be > 0")
    grid = dynamic_cast<grd_2d_xz_arakawa_c<real_t>*>(g);
    if (grid == NULL) error_macro("this version of the MPDATA scheme works with the Arakawa-C grid only!")
  }

  public: void op_1D(Array<quantity<unit, real_t>, 3>* psi[], 
    const Range &i, const Range &j, const Range &k, 
    const int n, const int step,
    const Array<quantity<si::dimensionless, real_t>, 3> &Cx, 
    const Array<quantity<si::dimensionless, real_t>, 3> &, 
    const Array<quantity<si::dimensionless, real_t>, 3> &Cz
  )
  {
#    define F_donorcl(psi_l,psi_r,U) (.5 * (U + sqrt(U*U)) * psi_l + .5 * (U - sqrt(U*U)) * psi_r)
#    define V_antydif(psi_l,psi_r,U) ((sqrt(U*U) - pow(U,2)) * (psi_r - psi_l) / (psi_r + psi_l))
    switch (step)
    {
      case 1:
      { 
        (*psi[n+1])(i,j,k) -= (
          F_donorcl((*psi[n])(i  ,j,k), (*psi[n])(i+1,j,k), Cx(i + grid->p_half,j,k)) - 
          F_donorcl((*psi[n])(i-1,j,k), (*psi[n])(i  ,j,k), Cx(i - grid->m_half,j,k))
        );
        break;
      }
      default:
      {
        (*psi[n+1])(i,j,k) -= (
          F_donorcl((*psi[n])(i,j,k), (*psi[n])(i+1,j,k), 
            where(
              (*psi[n])(i+1,j,k) + (*psi[n])(i,j,k) > 0,
              V_antydif((*psi[n])(i,j,k), (*psi[n])(i+1,j,k), Cx(i + grid->p_half,j,k)),
              quantity<unit, real_t>(0.)
            )
          ) -
          F_donorcl((*psi[n])(i-1,j,k), (*psi[n])(i,j,k), 
            where(
              (*psi[n])(i,j,k) + (*psi[n])(i-1,j,k) > 0,
              V_antydif((*psi[n])(i-1,j,k), (*psi[n])(i,j,k), Cx(i - grid->m_half,j,k)),
              quantity<unit, real_t>(0.)
            )
          )
        );
      }
    }
#    undef F
#    undef C
  }
};
#endif
