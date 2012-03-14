/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date February - March 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef RHS_HPP
#  define RHS_HPP

#  include "cmn.hpp" 
#  include "mtx.hpp"

/** @brief TODO
 *
 * a right-hand-side (RHS) of a generalised transport equation: 
 *
 * \f$
 *   \partial_t \psi + \nabla \cdot (\vec{v} \psi) = RHS_1 + RHS_2 + \ldots
 * \f$
 *
 * may be decomposed into a linear and nonlinear part with respect to \f$\psi\f$:
 *
 * \f$
 *   RHS = (C_1 + C_2 + \ldots) \cdot \psi + (R_1 + R_2 + \cdots)
 * \f$
 *
 * Assuming \f$C_i\f$ are constant in time and \f$R_i^{n+1}\f$ may be evaluated,
 * an NFT numerical integration procedure for advancing \f$\psi\f$ from 
 * time-level \f$n\f$ to time-level \f$n+1\f$ (separated by a timestep \f$\Delta t\f$)
 * may be then expressed as:
 *
 * \f$
 *   \psi^{n+1} = Src\left(
 *     Adv(Src(\psi^n, \sum_i R_i^{n}, \sum_i C_i), V^{n+1/2}),
 *     \sum_i R_i^{n+1}, \sum_i C_i
 *   \right) 
 * \f$
 *
 * where Adv is an advection operator (upstream, MPDATA, etc) and Src is a 
 * function applying the forcing terms:
 *
 * \f$
 *   Src(\psi, R, C) = \left[
 *     \psi \cdot + \frac{\Delta t}{2} R
 *   \right] \cdot \left[
 *     1 - \frac{\Delta t}{2} C
 *   \right]^{-1}
 * \f$
 *
 * the procedure takes advantage of the fact that the linear terms may be 
 * solved implicitly.
 *
 * cf. the discussion of implicit (linear) and explicit terms of the rhs in
 * section 2.2 of Prusa, Smolarkiewicz & Wyszogrodzki 2008 (Computers & Fluids)
 * and references therein...
 */
template <typename real_t>
class rhs : root
{
  // TODO: stencil-extent-like method?

  public: virtual void explicit_part(
    mtx::arr<real_t> &R, 
    mtx::arr<real_t> **aux, // TODO: const
    mtx::arr<real_t> **psi // TODO const
  ) 
  {} // TODO: better pure virtual?

  public: virtual real_t implicit_part(
    quantity<si::time, real_t> dt // TODO: const
  )
  { // TODO: better pure virtual?
    return real_t(0.);
  }
};  
#endif
