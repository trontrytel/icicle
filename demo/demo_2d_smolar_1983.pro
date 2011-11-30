; @copyright University of Warsaw
;  @date November 2011
;  @section LICENSE
;    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
;  @section DESCRIPTION
;    example usage: 
;    GDL> demo_2D, nx=20, nz=30, nt=100, adv='mpdata', el='uniform --vel.uniform.u .2 --vel.uniform.v 0 --vel.uniform.w .3'
;    GDL> demo_2D, nx=20, nz=30, nt=100, adv='mpdata', vel='rasinski --vel.rasinski.A 1 --vel.rasinski.Z_top 30 --vel.rasinski.Z_clb 20'
;    GDL> demo_2D, nx=100, nz=100, nt=628, adv='mpdata', vel='test --vel.test.omega .1', ini='cone --ini.cone.h 3.87 --ini.cone.x0 75 --ini.cone.z0 50 --ini.cone.r 15'
;    GDL> demo_2D, nx=50, nz=50, nt=200, adv='mpdata', vel='test --vel.test.omega .1', ini='cone --ini.cone.h 3.87 --ini.cone.x0 40 --ini.cone.z0 20 --ini.cone.r 5'

pro demo_2d_smolar_1983, test=test
  nt=6*628
  nx=100
  nz=100
  file='demo_2d_smolar_1983.nc'
  spawn, 'rm -f '+file
  cmd = '../icicle' + $
    ' --bits 32 --dt .1 --grd.dx 1 --grd.dy 1 --grd.dz 1'+ $
    '  --grd arakawa-c --ny 1' + $
    ' --nt '+string(nt)+ $
    ' --nx '+string(nx)+ $
    ' --nz '+string(nz)+ $
    ' --vel test --vel.test.omega .1'+ $
    ' --adv mpdata'+ $
    ' --ini cone --ini.cone.h 3.87 --ini.cone.x0 75 --ini.cone.z0 50 --ini.cone.r 15'+ $
    ' --slv serial' + $
    ' --out netcdf --out.netcdf.file '+file
  spawn, cmd
  a = ncdf_open(file) 
  scl = 2
  window, xsize=scl*nx, ysize=scl*nz
  loadct, 1
  i=0l
  while i lt nt-1 do begin ; repeat 100 times...
    ncdf_varget, a, 'psi', p, offset=[0,0,0,i++ ], count=[nz, 1, nx, 1]
    p = reform(p[*,0,*,0]) ; removing dimensions of size one
    p = rotate(p, 4)
;    set_plot, 'svg'
;    device, filename=string(i, FORMAT='(I05)')+".svg"
;    surface, p
;    device, /close
;stop
    tvscl, rebin(p, scl*nx, scl*nz, /sample)
  endwhile
end
