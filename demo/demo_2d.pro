;; @file
;  @author Sylwester Arabas (University of Warsaw) <slayoo@igf.fuw.edu.pl>
;  @date November 2011
;  @section LICENSE
;    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
;  @section DESCRIPTION
;    example usage: 
;    GDL> demo_2D, nx=20, nz=30, u=.2, w=.3, nt=100, adv='upstream'

pro demo_2D, nx=nx, nz=nz, u=u, w=w, nt=nt, adv=adv
  spawn, 'rm -f demo_2d.nc'
  spawn, '../icicle' + $
    ' --bits 32 --dt 1 --dx 1 --dy 1 --dz 1 --v 0 --ny 1' + $
    ' --nt ' +  strmid(nt, 2) + $
    ' --nx ' +  strmid(nx, 2) + $
    ' --nz ' +  strmid(nz, 2) + $
    ' --u '  +  strmid(u,  2) + $
    ' --w '  +  strmid(w,  2) + $
    ' --adv ' + adv + $
    ' --dom serial --nsd 2' + $
    ' --out netcdf --outfile demo_2d.nc'
  a = ncdf_open('demo_2d.nc') 
  scl = 20
  window, xsize=scl*nx, ysize=scl*nz
  loadct, 1
  i=0l
  while i lt 100 * nt do begin 
    ncdf_varget, a, 'psi', p, offset=[0,0,0,i++ mod nt], count=[nz, 1, nx, 1]
    p = reform(p[*,0,*,0]) ; removing dimensions of size one
    p = rotate(p, 4)
    tvscl, rebin(p, scl*nx, scl*nz, /sample)
  endwhile
end
