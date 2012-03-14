## @file
#  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
#  @copyright University of Warsaw
#  @date December 2011
#  @section LICENSE
#    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)

import numpy as np                       # arrays
from scipy.io.netcdf import netcdf_file  # netcdf IO
import subprocess                        # shell calls
import os                                # unlink()
from enthought.mayavi import mlab        

mlab.figure()
for adv in ('mpdata','upstream','leapfrog'):
  print 'running the model...'
  file = 'test_mayavi_'+adv+'.nc'
  freq = 1
  nx = 24#12
  ny = 16#8
  nz = 12#6
  nt = 50#50
  Cx = -.3
  Cy = -.5
  Cz = -.1
  cmd = (
    '../../icicle',
    '--dt','1', 
    '--grd.dx','1',
    '--grd.dy','1',
    '--grd.dz','1',
    '--nt',str(nt),
    '--grd.nx',str(nx),
    '--grd.ny',str(ny),
    '--grd.nz',str(nz),
    '--vel','uniform',
    '--vel.uniform.u',str(Cx),
    '--vel.uniform.v',str(Cy),
    '--vel.uniform.w',str(Cz),
    '--adv',adv,'--adv.mpdata.fct','0','--adv.mpdata.iord','4',
    '--ini','boxcar', 
    '--ini.boxcar.ax',str(nx-8), 
    '--ini.boxcar.ay',str(ny-6), 
    '--ini.boxcar.az',str(nz-4), 
    '--ini.boxcar.bx',str(nx), 
    '--ini.boxcar.by',str(ny), 
    '--ini.boxcar.bz',str(nz), 
    '--ini.boxcar.A','2', 
    '--ini.boxcar.A0','1', 
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
  x, y, z = np.mgrid[.5:nx:1,.5:ny:1,.5:nz:1]

  print 'saving animation frames...'
  for t in range(psi.shape[0]) :
    data = psi[t,:,:,:].copy()
    print data.min(), data.max()
    mlab.points3d(x, y, z, data.byteswap(), vmin=.6, vmax=3.4, colormap='spectral', scale_factor=.15, opacity=.9, reset_zoom=False)
    mlab.colorbar(orientation='horizontal')
    mlab.outline(extent=[0,nx,0,ny,0,nz])
    mlab.axes()
    #mlab.orientation_axes()
    mlab.text(.05,.9,'Courant = ('+str(Cx)+','+str(Cy)+','+str(Cz)+')      '+adv+" t/dt="+format(t*freq,"04d"),width=.9)
    mlab.get_engine().scenes[0].scene.isometric_view()
    mlab.move(forward=-6.5,up=-3.5)
    mlab.savefig('fig1-'+adv+'-'+format(t*freq,"05d")+'.png')
    mlab.clf()

  print 'closing the file...'
  nc.close()

  print 'composing animation...'
  cmd = ('convert','*fig1-'+adv+'-*.png','fig-'+adv+'-3d.gif')
  subprocess.check_call(cmd)

  print 'deleting frames and the netcdf file...'
  # TODO
