#!/bin/bash
rm tmp.nc; 
../../icicle --dt 1 --grd.dx 1 --grd.dy 1 --grd.dz 1 --nt 1 --bits 32 \
   --nx 3 --ny 6 --nz 4 \
  --vel uniform --vel.uniform.u -0.3 --vel.uniform.v -0.5 --vel.uniform.w -0.1 \
  --adv mpdata --adv.mpdata.fct 1 --adv.mpdata.iord 4 --adv.mpdata.cross_terms 1 \
  --ini boxcar \
    --ini.boxcar.ax 1 --ini.boxcar.ay 0 --ini.boxcar.az 1 \
    --ini.boxcar.bx 3 --ini.boxcar.by 5 --ini.boxcar.bz 2 \
    --ini.boxcar.A 1.5 --ini.boxcar.A0 1.0 \
  --slv serial \
  --out netcdf --out.netcdf.ver 3 --out.netcdf.freq 1 --out.netcdf.file tmp.nc 
ncdump tmp.nc | fgrep 0.999
