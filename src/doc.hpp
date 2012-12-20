/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date February 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 @  @brief contains the main page of the documentation
 */
/** @mainpage
 *  @section sec_ABOUT_ICICLE About icicle
 *
 * Icicle is an object-oriented C++ implementation of a nonoscillatory forward in time (NFT)
 * solver for systems of generalised transport equations with emphasis on cloud 
 * and precipitation modelling applications. 
 * Icicle solves systems of general transport equations of the following type:
 *
 * \f$ \partial_t \psi + \nabla \cdot (\vec{v} \psi) = R \f$
 * 
 * where
 *
 * \f$ \psi = [\psi_1, \psi_2, \ldots ] \f$ is a set of conservative dependent variables, 
 * \f$ R = [R_1, R_2, \ldots ] \f$ are the forcing terms, 
 * and \f$ \vec{v} = [u, v, w] \f$ is the velocity field.
 *
 * In most numerical aspects icicle follows the design of the MPDATA-based NFT
 * solvers of Smolarkiewicz et al., for a review and list of references consult:
 *           
 * Smolarkiewicz, P. and Margolin, L., 1998:
 * MPDATA: A Finite-Difference Solver for Geophysical Flows.
 * Journal of Computational Physics, 140, 459-480 
 *
 * For a list of required and optional packages and installation instructions 
 * consult the @ref sec_README
 *
 * Icicle's source code is publicly available under the terms of the @ref sec_COPYING
 * from a git repository at Github: http://github.com/slayoo/icicle/
 *
 * @section sec_EQUATIONS Supported equation sets
 * @table
 * | definition              | equation set | example |
 * ----------------------------------------------------
 * | eqs_scalar_advection    |              |         |
 * | eqs_harmonic_oscillator |              |         |
 * @endtable
 *
 * @section sec_CREDITS Credits
 *         
 * Icicle is an academic project developed by Sylwester Arabas, Anna Jaruga, 
 * and Hanna Pawlowska at the 
 * [Institute of Geophysics](http://www.igf.fuw.edu.pl/),
 * [Faculty of Physics](http://www.fuw.edu.pl/),
 * [University of Warsaw](http://www.uw.edu.pl/) (the copyright holder)
 * with funding from the Polish [National Science Centre](http://www.ncn.gov.pl/)
 * (cf. @ref sec_README for details).
 */
/** @page README README file (incl. requirements and installation instructions)
 *  @section sec_README README file
 *  @verbinclude "../README"
 */
/** @page HACKING HACKING file (coding conventions)
 *  @section sec_HACKING HACKING file (coding conventions)
 *  @verbinclude "../HACKING"
 */
/** @page COPYING GNU General Public License version 3
 *  @section sec_COPYING GNU General Public License version 3
 *  @verbinclude "../COPYING"
 */
