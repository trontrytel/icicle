## @file
#  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
#  @copyright University of Warsaw
#  @date Januar 2012
#  @section LICENSE
#    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)

import numpy as np                       # arrays
from scipy.io.netcdf import netcdf_file  # netcdf IO
import subprocess                        # shell calls
import os                                # unlink()
import sys                               # exit

file = 'tmp.nc'
for iord in (2,3,4): 
  print 'running 3D MPDATA-FCT simulation with iord=',str(iord),' ...'
  min = 1.
  max = 3.
  digits = 5

  freq = 1
  nx = 12
  ny = 8
  nz = 6
  nt = 5
  Cx = -.2
  Cy = -.3
  Cz = -.1
  dt = 1
  cmd = (
    '../../icicle','--bits','32',
    '--dt',str(dt), 
    '--grd.dx','1',
    '--grd.dy','1',
    '--grd.dz','1',
    '--t_max',str(nt * dt),
    '--grd.nx',str(nx),
    '--grd.ny',str(ny),
    '--grd.nz',str(nz),
    '--vel','uniform',
    '--vel.uniform.u',str(Cx),
    '--vel.uniform.v',str(Cy),
    '--vel.uniform.w',str(Cz),
    '--adv','mpdata',
      '--adv.mpdata.fct','1',
      '--adv.mpdata.cross_terms','1', #TEMP!
      '--adv.mpdata.third_order','0', 
      '--adv.mpdata.iord',str(iord),
    '--ini','boxcar', 
    '--ini.boxcar.ax',str(nx-8), 
    '--ini.boxcar.ay',str(ny-6), 
    '--ini.boxcar.az',str(nz-4), 
    '--ini.boxcar.bx',str(nx), 
    '--ini.boxcar.by',str(ny), 
    '--ini.boxcar.bz',str(nz), 
    '--ini.boxcar.A',str(max - min), 
    '--ini.boxcar.A0',str(min), 
    '--slv','threads','--nsd','2', 
    '--out','netcdf', 
    '--out.netcdf.ver','3', 
    '--dt_out',str(freq*dt),
    '--out.netcdf.file',file
  )
  if os.path.exists(file) : os.unlink(file)
  subprocess.check_call(cmd)

  nc = netcdf_file(file)
  psi = nc.variables['psi']

  for t in range(psi.shape[0]) :
    data = psi[t,:,:,:]
    print "data.min()-min:", data.min()-min, "data.max()-max", data.max()- max
    if (round(data.min(),digits) < min) or (round(data.max(),digits) > max) :
      print data.min(), " < ", min, " or ", data.max(), " > ", max
      sys.exit(1) 

  print 'closing and deleting the file...'
  nc.close()
  if os.path.exists(file) : os.unlink(file)
