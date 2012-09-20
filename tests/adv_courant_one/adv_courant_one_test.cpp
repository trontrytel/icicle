/** @file
 *  @example adv_courant_one/adv_courant_one_test.cpp
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011 - September 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @section DESCRIPTION
 *    sanity checks for one-dimensional advection with Courant
 *    number of unity (should result in exact solution)
 */

#include <cstdlib>

#include <list>
using std::list;

#include <sstream>
using std::string;
using std::ostringstream;

#include <iostream>
using std::cerr;
using std::endl;

#include "../../src/cfg/cfg_types.hpp"
#include "../../src/cfg/cfg_boost_thread.hpp"

int main()
{
  for (int &bits : list<int>({
#if defined(USE_FLOAT)
    8 * sizeof(float),
#endif
#if defined(USE_DOUBLE)
    8 * sizeof(double),
#endif
#if defined(USE_LDOUBLE)
    8 * sizeof(long double),
#endif
#if defined(USE_FLOAT128)
    8 * sizeof(__float128),
#endif
  })) 
  {
    for (int &u : list<int>({-1,1}))
    {
      for (int &t_max : list<int>({10,30})) 
      {
        for (string &slv : list<string>({
          "serial",
#if defined(_OPENMP)
          "openmp",
#endif
#if defined(USE_BOOST_THREAD)
          "threads",
#endif
        })) // TODO: fork+openmp, fork (output not ready), fork+threads, mpi
        {
          for (int &nsd : list<int>({1,2,4,5}))
          {
            for (string &nxyz : list<string>({
              "--grd.nx 20 --vel.uniform.u",
              "--grd.ny 20 --vel.uniform.v",
              "--grd.nz 20 --vel.uniform.w",
            }))
            {
              for (string &adv : list<string>({
                "mpdata",
                "mpdata --adv.mpdata.iord 1",
                "mpdata --adv.mpdata.iord 3",
                "leapfrog",
                "mpdata --adv.mpdata.fct 1",
              }))
              {
                // combinations that does not make sense
                if (slv == "serial" and nsd == 1) continue; 
                if (nxyz != "--grd.nx 20 --vel.uniform.u" && nsd > 1) continue; // parallelism only in X
/*
// TODO!
                  if [ "$adv" = "mpdata --adv.mpdata.fct 1" -a $nsd = 20 ]; then continue; fi # halo > nx

*/
                ostringstream cmd;
                cmd << "../../icicle --ini boxcar --ini.boxcar.bx 1";
                cmd << " --bits " << bits << " --dt_out 1 --grd.dx 1 --grd.dy 1 --grd.dz 1";
                cmd << " --vel uniform " << nxyz << " " << u;
                cmd << " --t_max " << t_max << " --adv " << adv << " --slv " << slv << " --out gnuplot --nsd " << nsd;
                cmd << " 2>/dev/null"; 
                cmd << " | awk 'BEGIN {if (NF==0) print \"error\"; RS=\"   \"} {print}' | tail -21 | tr \"\\n\" \" \" ";
                cmd << " | awk '{if ($0!=\"0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0  \") exit 1;}'";

                int status = system(cmd.str().c_str());
                if (status != 0) 
                {
                  cerr << "KO: " << cmd.str() << endl;
                  exit(EXIT_FAILURE);
                }
                cerr << "OK: " << cmd.str() << endl;
              } // adv
            } // nxyz
          } // nsd
        } // slv
      } // t_max
    } // u
  } // bits
} // main()
