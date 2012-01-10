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
for iord in (1,2,3,4):
  print 'running 3D MPDATA-FCT simulation with iord=',str(iord),' ...'
  min = 1.
  max = 3.
  digits = 5

  freq = 1
  nx = 24#12
  ny = 16#8
  nz = 12#6
  nt = 5
  Cx = -.3
  Cy = -.5
  Cz = -.1
  cmd = (
    '../../icicle','--bits','32',
    '--dt','1', 
    '--grd.dx','1',
    '--grd.dy','1',
    '--grd.dz','1',
    '--nt',str(nt),
    '--nx',str(nx),
    '--ny',str(ny),
    '--nz',str(nz),
    '--vel','uniform',
    '--vel.uniform.u',str(Cx),
    '--vel.uniform.v',str(Cy),
    '--vel.uniform.w',str(Cz),
    '--adv','mpdata','--adv.mpdata.fct','1','--adv.mpdata.iord',str(iord),
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
    '--out.netcdf.freq',str(freq),
    '--out.netcdf.file',file
  )
  if os.path.exists(file) : os.unlink(file)
  subprocess.check_call(cmd)

  print 'opening the output file...' # TODO: move it into the loop!
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
