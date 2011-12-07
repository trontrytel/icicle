## @file
#  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
#  @copyright University of Warsaw
#  @date November 2011
#  @section LICENSE
#    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)

for bits in {32,64,128}; do
  for u in {-1,1}; do
    for nt in {10,30}; do
      for slv in `../../icicle --slv list`; do
#        for sdom in `../../icicle --sdom list`; do
          for nsd in {1,2,4,5,20}; do
#            for nssdom in {1,2,3}; do
              for nout in {1,5,10}; do
                for adv in {mpdata,"mpdata --adv.mpdata.iord 1","mpdata --adv.mpdata.iord 3",leapfrog,"mpdata --adv.mpdata.fct 1"}; do
                  # combinations that does not make sense
                  if [ $slv = "serial" -a $nsd != 1 ]; then continue; fi
                  if [ $slv = "fork" ]; then continue; fi # TODO: output not ready
                  if [ $slv = "fork+openmp" ]; then continue; fi # TODO: not ready
                  if [ $slv = "fork+threads" ]; then continue; fi # TODO: not ready
#                  if [ $sdom = "serial" -a $nssdom != 1 ]; then continue; fi

                  # the actual test
                  cmd="../../icicle --ini boxcar --ini.boxcar.b 1"
                  cmd="$cmd --bits $bits --dt 1 --grd.dx 1 --grd.dy 1 --grd.dz 1 --grd arakawa-c-lorenz"
                  cmd="$cmd --vel uniform --vel.uniform.u $u --vel.uniform.v 0 --vel.uniform.w 0"
                  cmd="$cmd --nt $nt --nx 20 --ny 1 --nz 1 --adv $adv --slv $slv --out gnuplot --nsd $nsd"
# TODO: --sdom $sdom --nssdom $nssdom --nout $nout"
                  if [ $slv = "mpi" ]; then cmd="openmpirun -np $nsd $cmd"; fi;
                  ret=`$cmd 2>/dev/null` 
                  if [ $? != "0" ]; then 
                    echo "exited with non-zero status for: " $cmd
                    exit 1; 
                  fi 
                  ret=`echo "$ret" | awk 'BEGIN {RS="   "} {print}' | tail -21 | tr "\n" " "`
                  if [ "$ret" != '0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0  ' ]; then 
                    echo "wrong result: $ret" 
                    echo "command: " $cmd
                    exit 1; 
                  fi
                  echo "OK for: " $cmd
                done # adv
              done # nout
#            done # nssdom
          done # nsd
#        done # sdom
      done # slv
    done # nt
  done # u
done # bits