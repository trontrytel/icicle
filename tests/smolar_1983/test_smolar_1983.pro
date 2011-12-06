;  @copyright University of Warsaw
;  @date November 2011
;  @section LICENSE
;    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)

pro test_smolar_1983, full=full
  nt=80 ;6*628 ; TEMP!!!
  nx=100
  nz=100
  file='test_smolar_1983.nc'
;  spawn, 'rm -f '+file
  cmd = '../../icicles' + $
    ' --dt .1 --grd.dx 1 --grd.dy 1 --grd.dz 1'+ $
    ' --nt '+string(nt)+ $
    ' --nx '+string(nx)+ $
    ' --nz '+string(nz)+ $
    ' --vel test --vel.test.omega .1'+ $
    ' --adv mpdata --adv.mpdata.iord 1'+ $
    ' --ini cone --ini.cone.h 3.87 --ini.cone.x0 75 --ini.cone.z0 50 --ini.cone.r 15'+ $
    ' --slv threads --nsd 2' + $
    ' --out netcdf --out.netcdf.ver 3 --out.netcdf.freq 1 --out.netcdf.file '+file
;  spawn, cmd, exit_status=status
;  if status ne 0 then exit, status=1
  f = ncdf_open(file) 
  ncdf_varget, f, 'psi', psi
  for i=0,nt -1 do begin
    surface, reform(psi[*,0,*,i]), zrange=[0,4], title='MPDATA t=' + strtrim(nt, 2)
  print, i, total(psi[*,0,*,i])
  endfor
  ncdf_close, f
end
