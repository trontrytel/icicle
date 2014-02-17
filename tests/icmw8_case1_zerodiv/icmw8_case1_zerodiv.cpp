#include "../../src/icmw8_case1.hpp"

#include <libmpdata++/solvers/mpdata.hpp>
#include <libmpdata++/concurr/serial.hpp>
using namespace libmpdataxx;

int main()
{
  struct ct_params_t : ct_params_default_t
  {
    using real_t = float;
    enum { n_dims = 2 }; 
    enum { n_eqs = 2 };
    enum { opts = formulae::opts::nug };
    struct ix { enum {th, rv}; };
  };

  using slv_t = solvers::mpdata<ct_params_t>;

  slv_t::rt_params_t p;
  int nx=11, ny=11;
  p.span = {nx, ny};
 
  concurr::serial<slv_t, bcond::cyclic, bcond::open> run(p);
  
  icmw8_case1::intcond(run);

  blitz::Array<typename ct_params_t::real_t, 2> div(nx, ny); 

  blitz::Range i(0, nx-1), j(0, ny-1);
  div = (
    run.advector(0)(i,   j  ) - 
    run.advector(0)(i+1, j  )
  ) + ( 
    run.advector(1)(i,   j  ) - 
    run.advector(1)(i,   j+1)
  );

  std::cerr << run.advector(0) << std::endl;
  std::cerr << std::endl;
  std::cerr << run.advector(1) << std::endl;
  std::cerr << std::endl;
  std::cerr << div << std::endl;
  std::cerr << "eps = " << blitz::epsilon(typename ct_params_t::real_t(0)) << std::endl;

}
