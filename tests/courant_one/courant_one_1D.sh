## @file
#  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
#  @copyright University of Warsaw
#  @date November 2011
#  @section LICENSE
#    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)

for bits in {32,64,128}; do
  for u in {-1,1}; do
    for t_max in {10,30}; do
      for slv in `../../icicle --slv list`; do
#        for sdom in `../../icicle --sdom list`; do
          for nsd in {1,2,4,5}; do
            for nxyz in {"--grd.nx 20 --vel.uniform.u","--grd.ny 20 --vel.uniform.v","--grd.nz 20 --vel.uniform.w"}; do
#            for nssdom in {1,2,3}; do
              for nout in {1,5,10}; do
                for adv in {mpdata,"mpdata --adv.mpdata.iord 1","mpdata --adv.mpdata.iord 3",leapfrog,"mpdata --adv.mpdata.fct 1"}; do
                  # combinations that does not make sense
                  if [ $slv = "serial" -a $nsd != 1 ]; then continue; fi
                  if [ $slv = "fork" ]; then continue; fi # TODO: output not ready
                  if [ $slv = "fork+openmp" ]; then continue; fi # TODO: not ready
                  if [ $slv = "fork+threads" ]; then continue; fi # TODO: not ready
                  if [ "x$nxyz" \> "x--grd.nxx" -a $nsd -gt 1 ]; then continue; fi # parallelism only in X
                  if [ "$adv" = "mpdata --adv.mpdata.fct 1" -a $nsd = 20 ]; then continue; fi # halo > nx

                  # the actual test
                  cmd="../../icicle --ini boxcar --ini.boxcar.bx 1"
                  cmd="$cmd --bits $bits --dt_out 1 --grd.dx 1 --grd.dy 1 --grd.dz 1 --grd arakawa-c-lorenz"
                  cmd="$cmd --vel uniform $nxyz $u"
                  cmd="$cmd --t_max $t_max --adv $adv --slv $slv --out gnuplot --nsd $nsd"
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
            done # nxyz
          done # nsd
#        done # sdom
      done # slv
    done # nt
  done # u
done # bits
