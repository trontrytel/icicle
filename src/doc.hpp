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
 *           Icicle is a modern C++ implementation of a nonoscillatory forward in time (NFT)
 *           solver for systems of generalised transport equations with emphasis on cloud 
 *           and precipitation modelling applications. 
 *           It is being designed as a tool for integrating equations in one, two or three dimensions.
 *
 *           In most numerical aspects icicle follows the design of the MPDATA-based NFT
 *           solvers of Smolarkiewicz et al., for a review and list of references consult:
 *           
 *           Smolarkiewicz, P. and Margolin, L., 1998:<br />
 *           MPDATA: A Finite-Difference Solver for Geophysical Flows.<br />
 *           Journal of Computational Physics, 140, 459-480 
 *
 *           Icicle code is developed with the aim of providing a tool that:
 *           - is reliable, mainainable and auditable<br />
 *             (generation of reproducible results for research purposes),
 *           - has well documented, modular and human-readible code<br />
 *             (easy to study and extend for [PhD] students and researchers alike)
 *           - provides high performace and scalability in diverse parallel enviroments 
 *             without sacrificing code readibility or ease of use<br />
 *             (full separation of physics, numerics and technical aspects of the code),
 * 
 *           To achieve it icicle uses modern coding techniques and bases a lot of its functionality 
 *           on other free-libre-open-source software:
 *           - Blitz++ library is used to achieve better-than-FORTRAN performance with array-oriented number crunching
 *           - Boost.units library is used to enable dimensional analysis of the program code at compile-time
 *           - Boost.Thread and Boost.MPI are used for parallelisation
 *           - Boost.program_options, Boost.conversion, Boost.timer and others are used for what their names suggest
 *           - NetCDF library is used for both input and output
 *           - CMake and Doxygen are used for pre-compilation configuration and documentation-generation, respectively
 *           - last but not list, the GNU OpenMP-enabled C++11 compiler and other toold from the GNU developement chain are used
 *
 *           For a full list of required and optional packages and installation instructions 
 *           consult the @ref sec_README
 *
 *           Icicle's source code is publicly available under the terms of the @ref sec_COPYING
 *           from a git repository at Github: http://github.com/slayoo/icicle/
 *
 *  @section sec_CREDITS Credits
 *         
 *           Icicle is an academic project developed by Sylwester Arabas, Anna Jaruga, 
 *           and Hanna Pawlowska at the 
 *           <a href="http://www.igf.fuw.edu.pl/">Institute of Geophysics</a>,
 *           <a href="http://www.fuw.edu.pl/">Faculty of Physics</a>, 
 *           <a href="http://www.uw.edu.pl/">University of Warsaw</a> (copyright holder)
 *           with funding from the Polish <a href="http://www.ncn.gov.pl/">National Science Centre</a> 
 *           (cf. @ref sec_README for details).
 *
 *  @section sec_ABOUT_DOCS About this documentation
 *           
 *           This documentation is intended to serve both a Users' and a Developer's manual.
 *           It also covers presentation of the simulation results obtained by running the
 *           icicle's test suite. The test suite programs are listed in the <a href="examples.html">Examples</a>
 *           are intended to serve as a tutorial as well. 
 *
 *  @section sec_PROBLEM The mathematical/numerical problem and selected corresponding program components
 * 
 *           Icicle is being developed to serve as a numerical solver for systems of general transport
 *           equations of the following type:
 *
 *           \f$ \partial_t \psi + \nabla \cdot (\vec{v} \psi) = R \f$
 * 
 *           where
 *
 *           \f$ \psi = [\psi_1, \psi_2, \ldots ] \f$, \f$ R = [R_1, R_2, \ldots ] \f$ 
 *           and \f$ \vec{v} = [u, v, w] \f$
 *
 *           Icicle is written in an object-oriented manner and hence numerous
 *           terms and entities from the mathematical formulation of the model
 *           heve their counterparts in the program's class hierarchy.
 *           In the case of the above transport equation:
 *           - the @ref rhs class represents the \f$R\f$ term (right-hand-side, source term)
 *           - the @ref adv class represents the logic to numerically represent the advective 
 *             term \f$\nabla \cdot (\vec{v}\psi)\f$
 *           - the @ref slv class represents the time-stepping procedures used to advence the 
 *             solution in time (i.e. integration corresponding to the \f$\partial_t\f$ derivative)
 *           - the @ref vel class represents the velocity field \f$\vec{v}\f$
 *
 *           All of the above classes are pure vitual in C++ nomenclature, meaning that they
 *           do not contain particular algorithms or data structures but rather the requirements
 *           needed to be fullfilled in order to serve as a valid instance (object) of that class.
 *           For example, there are several numerical representations of the advection operator available:
 *           - the upstream algorithm represented by the adv_upstream class (inheriting from @ref adv)
 *           - the MPDATA algorithm represented by the adv_mpdata class (inheriting from adv_upstream)
 *           - the MPDATA-FCT algorithm represented by the adv_mpdata_fct class (inheriting from adv_mpdata)
 *           - the leapfrog algorithm represented by the adv_leapfrog class (inheriting from @ref adv)
 * 
 *  @section sec_OPTIONS The command-line interface 
 *
 *           All program parameters (including the choice of advection algorithm from the list above) 
 *           are controlled with UNIX command-line options. A typical calling sequence looks like:
 *           
 *           <tt>
 *           \$ icicle --bits 64 \\                                    <br /> 
 *             --eqs shallow_water \\                                  <br />
 *             --adv mpdata --adv.mpdata.iord 2 \\                     <br />
 *             --ini netcdf --ini.netcdf.file ini.nc \\                <br />
 *             --grd.dx 1 --grd.nx 100 \\                              <br />
 *             --nt 1000 --dt 1 --nout 10 \\                           <br />
 *             --out netcdf --out.netcdf.file out.nc \\                <br />
 *             --slv threads --nsd 4 \\                                <br />
 *             ...    
 *           </tt>
 *           
 *           The example command above would instruct the solver to integrate the 
 *           shallow-water equations (<tt>--eqs shallow_water</tt>) using the MPDATA 
 *           iterative algorithm (<tt>--adv mpdata</tt>) with three iterations
 *           (<tt>--adv.mpdata.iord 2</tt>). The computational grid would consist of
 *           100 points (<tt>--grd.nx 100</tt>) spaced every 1 meter (<tt>--grd.dx 1</tt>).
 *           The solver would use timestep of 1 second (<tt>--dt 1</tt>) and perform
 *           1000 integration steps (<tt>--nt 1000</tt>) outputting the results
 *           every 10-th step (<tt>--nout 10</tt>) to a netCDF file (<tt>--out netcdf</tt>)
 *           named out.nc (<tt>--out.netcdf.file out.nc</tt>). The computations would
 *           be performed in double precision (<tt>--bits 64</tt>), in parallel
 *           by using four (<tt>--nsd 4</tt>) memory-sharing threads (<tt>--slv threads</tt>).
 *
 *           A complete list of options may be obtained by executing <tt>\$ icicle --help</tt>.
 *
 *  @section sec_CONFLOW An outline of the control flow and program structure
 *
 *           Program execution starts in the main() function of the icicle.cpp file.
 *           There the command-line parsing begins by first verifying the correctness of
 *           the option names and then by reading the value of the <tt>--bits</tt> option
 *           and executing mdl<real_t>() function with the template parameter set to
 *           real_t=float for 32 bits, real_t=double for 64 bits etc, consequently
 *           enabling one to switch floating point precision without recompilation.
 *           
 *           The role of the mdl<real_t>() function is to:
 *           - trigger parsing of the remaining command-line options (by calling opt_adv(), opt_vel(), opt_ini(), opt_eqs(), ...), 
 *           - create a data structure representing the simullation set-up (an instance of the @ref stp structure) containing:
 *             - an instance of a @ref grd subclass - representing the chosen grid type (the default is grd_arakawa_c_lorenz),
 *             - an instance of a @ref adv subclass - representing the chosen advection algorithm,
 *             - an instance of a @ref vel subclass - representing the chosen velocity fiels,
 *             - an instance of a @ref ini subclass - representing the chosen initial condition,
 *             - an instance of a @ref eqs class - representing the chosen system of equations,
 *             - and manually or automatically chosen timestep stp::dt 
 *           - instantiate an object responsible for storing the results (e.g. the out_netcdf),
 *           - instantiate a solver (an object of a @ref slv subclass)  initialising it with the set-up structure and associating it with the output and 
 *           - start the integration loop of the solver (the @ref slv::integ_loop() method).
 *
 *           The logic behind the integration loop is in the slv_parallel::integ_loop_sd() method
 *           (applicable to both serial and parallel computations, despite the name).
 *
 *  @section sec_FEATURES List of key implemented and planned features
 *           The feature list below is intentionally not exhaustive and does not contain 
 *           for example several modules intended for debugging only.
 *           Features not fully implemented yet or planned to be implemented soon are typed with italics.
 *
 *           - <b>solver operation modes:</b> (the <tt>--slv</tt> option)
 *             - serial (slv_parallel_serial)
 *             - parallel with shared-memory (slv_parallel)
 *               - OpenMP-based (slv_parallel_openmp)
 *               - Boost.Thread-based (slv_parallel_threads)
 *               - <i>C++11 threads-based</i>
 *             - parallel with non-shared memory (slv_parallel_distmem)
 *               - MPI-based (slv_parallel_distmem_mpi)
 *               - fork() + IPC communications (slv_prallel_distmem_fork)
 *             - <i>MPI+threads combinations</i>
 *           - <b>advection algorithms:</b> (the <tt>--adv</tt> option)
 *             - upstream (adv_upstream)
 *             - MPDATA with arbitrary number of iterations and optional 3-rd order accuracy (adv_mpdata)
 *             - Flux-Corrected (aka non-oscillatory) version of MPDATA (adv_mpdata_fct)
 *             - leapfrog (adv_leapfrog)
 *           - <b>grids:</b> (the <tt>--grd</tt> option):
 *             - Arakawa-C/Lorenz staggered grid (grd_arakawa-c-lorenz)
 *           - <b>initial condition specification</b> (the <tt>--ini</tt> option):
 *             - functional form - for single equation only (ini_func):
 *               - Gaussian shape (ini_func_gauss)
 *               - boxcar shape (ini_func_boxcar)
 *               - 2D cone (ini_func_cone, see also @ref test_smolar_1983.cpp)
 *             - data fed from a netCDF file - arbitrary number of variables (ini_netcdf)
 *           - <b>velocity field</b> (the <tt>--vel</tt> option):
 *             - functional form (kinematic simulations):
 *               - uniform (vel_func_uniform)
 *               - Stratocumulus-like eddy from Rasinski et al. 2011 (vel_func_rasinski) 
 *               - <i>8th WMO Cloud Modelling Workshop Case 1 velocity field</i>
 *             - momentum-equation derived:
 *               - extrapolated to n+1/2 using n and n-1 (vel_momeq_extrapol)
 *           - <b>equation systems</b> (the <tt>--eqs</tt> option):
 *             - homogeneous transport of a single scalar (eqs_scalar_advection)
 *             - shallow water equations (eqs_shallow_water)
 *             - <i>2D heat and moisture transport model</i> (8th WMO Cloud Modelling Workshop Case 1)
 *           - <i><b>moist thermodynamics related source terms</b></i>:
 *             - <i>Kessler parameterisation</i>
 *             - <i>Super-Droplet parameterisation</i>
 *           - <b>timestep choice:</b> (the <tt>--dt</tt> option, see constructors of the @ref stp structure):
 *             - automaticly chosen based on specified dt_out and the type of advection algorithm (for kinematic simulations only)
 *             - manually specified
 *           - <b>floating point precision choices</b> (the <tt>--bits</tt> option):
 *             - 32 bit aka "float"
 *             - 64 bit aka "double"
 *             - 80 bit (on most compilers on x86_64 architectures) aka "long double"
 *             - <i>128 bit using the __float128 type - a GCC extension</i>
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
