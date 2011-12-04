;  @copyright University of Warsaw
;  @date November 2011
;  @section LICENSE
;    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)

pro test_smolar_1983, full=full
  nt=10 ;6*628 ; TEMP!!!
  nx=100
  nz=100
  file='test_smolar_1983.nc'
  spawn, 'rm -f '+file
  cmd = '../../icicles' + $
    ' --dt .1 --grd.dx 1 --grd.dy 1 --grd.dz 1'+ $
    ' --nt '+string(nt)+ $
    ' --nx '+string(nx)+ $
    ' --nz '+string(nz)+ $
    ' --vel test --vel.test.omega .1'+ $
    ' --adv mpdata'+ $
    ' --ini cone --ini.cone.h 3.87 --ini.cone.x0 75 --ini.cone.z0 50 --ini.cone.r 15'+ $
    ' --slv serial' + $
    ' --out netcdf --out.netcdf.freq '+string(nt)+' --out.netcdf.file '+file
  spawn, cmd, exit_status=status
  if status ne 0 then exit, status=1
  f = ncdf_open(file) 
  ncdf_varget, f, 'psi', psi
  set_plot, 'svg'
  device, filename='fig1.svg'
  !P.MULTI=[0,2,1]
  surface, reform(psi[*,0,*,0]), zrange=[0,4], title='MPDATA (iord=2) t=0'
  surface, reform(psi[*,0,*,1]), zrange=[0,4], title='MPDATA (iord=2) t=' + strtrim(nt, 2)
  device, /close
  ncdf_close, f
end
