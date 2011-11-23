for bits in {32,64,128}; do
  for u in {-1,1}; do
    for nt in {10,30}; do
      for dom in `../icicle --dom list`; do
#        for sdom in `../icicle --sdom list`; do
          for nsd in {1,2,4,5,20}; do
            for nssdom in {1,2,3}; do
              for nout in {1,5,10}; do
                for adv in {mpdata,"mpdata --adv.mpdata.iord 1","mpdata --adv.mpdata.iord 3",leapfrog}; do
                  # combinations that does not make sense
                  if [ $dom = "serial" -a $nsd != 1 ]; then continue; fi
#                  if [ $dom = "serial" -a $sdom != "serial" ]; then continue; fi
#                  if [ $sdom = "serial" -a $nssdom != 1 ]; then continue; fi

                  # the actual test
                  cmd="../icicle"
                  cmd="$cmd --bits $bits --dt 1 --grd.dx 1 --grd.dy 1 --grd.dz 1 --grd arakawa-c"
                  cmd="$cmd --vel uniform --vel.uniform.u $u --vel.uniform.v 0 --vel.uniform.w 0"
                  cmd="$cmd --nt $nt --nx 20 --ny 1 --nz 1 --adv $adv --dom $dom --out gnuplot --nsd $nsd"
# TODO: --sdom $sdom --nssdom $nssdom --nout $nout"
                  if [ $dom = "mpi" ]; then cmd="openmpirun -np $nsd $cmd"; fi;
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
            done # nssdom
          done # nsd
#        done # sdom
      done # dom
    done # nt
  done # u
done # bits
