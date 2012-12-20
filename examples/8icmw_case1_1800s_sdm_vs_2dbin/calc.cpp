/* @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @author Ania Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date March 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#include <vector>
using std::vector;

#include <string>
using std::string;

#include <iostream>
using std::ostringstream;
using std::cerr;
using std::endl;

#include <boost/units/systems/si.hpp>
#include <boost/units/io.hpp>
namespace si = boost::units::si;
using boost::units::quantity;
using boost::units::divide_typeof_helper;
using boost::units::detail::get_value;

#include <boost/filesystem.hpp>

#include "../../src/cmn/cmn_error.hpp"
#include "../../src/cmn/cmn_units.hpp"
#include "../../src/phc/phc.hpp"
#include "../../src/wkc/wkc_icmw8_case1.hpp"

typedef double real_t;

#include "bins.hpp"

int main(int argc, char **argv)
{
  string dir = argc > 1 ? argv[1] : "tmp";
  notice_macro("creating " << dir << "...");
  boost::filesystem::create_directory(dir);

  notice_macro("creating input files ...")
  string opts = wkc::icmw8_case1::create_files(
    wkc::icmw8_case1::sdm,
    dir + "/rho.nc", 
    dir + "/ini.nc",
    50, // nx
    50, // ny
    real_t(30) * si::metres, // dx
    real_t(30) * si::metres  // dy
  );

  notice_macro("calling the solver ...")
  ostringstream cmd;
  cmd
    << "../../icicle " << opts 
    << " --adv mpdata"
      << " --adv.mpdata.fct " << false
      << " --adv.mpdata.iord " << 2
      << " --adv.mpdata.third_order " << false
    << " --dt " << real_t(1)
    << " --nt " << real_t(100) // 1800
    << " --nout " << real_t(5) // 60
    << " --out netcdf" 
    << " --out.netcdf.file " << dir << "/out.nc"
    << " --slv serial"

    << " --eqs.todo_sdm.adve.sstp " << 1
    << " --eqs.todo_sdm.sedi.sstp " << 1
    << " --eqs.todo_sdm.cond.sstp " << 5
    << " --eqs.todo_sdm.chem.sstp " << 0
    << " --eqs.todo_sdm.coal.sstp " << 10

    << " --eqs.todo_sdm.adve " << true
    << " --eqs.todo_sdm.cond " << true
    << " --eqs.todo_sdm.coal " << true
    << " --eqs.todo_sdm.sedi " << true
    << " --eqs.todo_sdm.chem " << false
    << " --eqs.todo_sdm.sd_conc_mean " << 128

    // Zach's output setting (mass-doubling layout)
    << " --eqs.todo_sdm.out_m0 \"";
  vector<quantity<si::length, real_t>> left_edges = bins();
  for (int i = 0; i < left_edges.size()-1; ++i) 
    cmd << real_t(left_edges[i] / si::metres) << ":" << real_t(left_edges[i + 1] / si::metres) << ";";
  cmd << "\"";
    
  if (EXIT_SUCCESS != system(cmd.str().c_str()))
    error_macro("model run failed: " << cmd.str())
}
