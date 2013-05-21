#include <libmpdata++/bcond/cyclic_2d.hpp>
#include <libmpdata++/concurr/threads.hpp>
#include <libmpdata++/output/gnuplot.hpp>

#include "kin_cloud_2d_blk_1m.hpp"

enum {x, z}; // dimensions
using real_t = double;

enum {rhod_th_ix, rhod_rv_ix, rhod_rc_ix, rhod_rr_ix }; // variables

// 8th ICMW case 1 by Wojciech Grabowski)
#include "icmw8_case1.hpp"

// simulation and output parameters
template <class T>
void setopts(T &params, int nt, int n_iters)
{
  params.outfreq = nt / 50; 
  //params.gnuplot_zrange = p.gnuplot_cbrange = "[.5:2.5]"; // TODO: per variable!
  params.gnuplot_view = "map";
  {
    std::ostringstream tmp;
    tmp << "output/figure_iters=" << n_iters << "_%s_%d.svg";
    params.gnuplot_output = tmp.str();
  }
  params.outvars = 
  {
    {rhod_th_ix, {.name = "\\rho_d \\theta", .unit = "kg/m^{-3} K"}},
    {rhod_rv_ix, {.name = "\\rho_v", .unit = "kg/m^{-3}"}},
    {rhod_rc_ix, {.name = "\\rho_c", .unit = "kg/m^{-3}"}},
    {rhod_rr_ix, {.name = "\\rho_r", .unit = "kg/m^{-3}"}}
  };

  // Kessler scheme options
  params.bulk_opts.cevp = true;
  params.bulk_opts.revp = true;
  params.bulk_opts.conv = true;
  params.bulk_opts.clct = true;
  params.bulk_opts.sedi = true;
}


int main()
{
  int nx = 32, nz = 32, nt = 500;
  const int n_iters = 2;

  // helper type to shorten the code below
  using solver_t = output::gnuplot<kin_cloud_2d_blk_1m<real_t, n_iters, solvers::strang,
    rhod_th_ix, 
    rhod_rv_ix,
    rhod_rc_ix,
    rhod_rr_ix
  >>;

  // instantiation of structure containing simulation parameters
  solver_t::params_t p;

  // output and simulation parameters
  setopts(p, nt, n_iters);
  icmw8_case1::setopts(p, nz);

  // solver instantiation
  concurr::threads<solver_t, bcond::cyclic, bcond::cyclic> slv(nx, nz, p);

  // initial condition
  icmw8_case1::intcond(slv);

  // timestepping
  slv.advance(nt);
}
