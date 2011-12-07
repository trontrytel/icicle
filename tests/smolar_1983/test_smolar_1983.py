import numpy as np                       # arrays
from scipy.io.netcdf import netcdf_file  # netcdf IO
import subprocess                        # shell calls
import os                                # unlink()
import matplotlib                        # plotting
matplotlib.use('Agg')                    # non-interactive backend
import matplotlib.pyplot as plt          # ?
from mpl_toolkits.mplot3d import Axes3D  # 3D plot
from matplotlib import cm                # colormap 

for adv in ('leapfrog','upstream','mpdata'):
  print 'running the model...'
  file = 'test_smolar_1983'+adv+'.nc'
  if os.path.exists(file) : os.unlink(file)
  freq = 25
  cmd = (
    '../../icicle',
    '--dt','.1', 
    '--grd.dx','1',
    '--grd.dy','1',
    '--grd.dz','1',
    '--nt','1256', # should be 6*628 !
    '--nx','100',
    '--nz','100',
    '--vel','test',
    '--vel.test.omega','.1',
    '--adv',adv,
    '--ini','cone', 
    '--ini.cone.h','3.87', 
    '--ini.cone.x0','75', 
    '--ini.cone.z0','50', 
    '--ini.cone.r','15',
    '--slv','threads','--nsd','2',
    '--out','netcdf', 
    '--out.netcdf.ver','3', 
    '--out.netcdf.freq',str(freq),
    '--out.netcdf.file',file
  )
  subprocess.check_call(cmd)

  print 'opening the output file...' # TODO: move it into the loop!
  nc = netcdf_file(file)
  psi = nc.variables['psi']

  print 'setting-up the plot...'
  x = np.arange(.5, 100.5, 1)
  y = np.arange(.5, 100.5, 1)
  x, y = np.meshgrid(x, y)

  print 'saving animation frames...'
  fig = plt.figure()
  vmin = -.5
  vmax = 4.
  for t in range(psi.shape[0]) :
    axs = fig.gca(projection='3d')
    srf = axs.plot_surface(x, y, psi[t,:,0,:], rstride=1, cstride=1, cmap=cm.gist_ncar, linewidth=0, antialiased=True, vmin=vmin, vmax=vmax)
    axs.set_zlim3d(vmin, vmax)
    fig.colorbar(srf) #shrink=0.5, aspect=5)
    fig.suptitle(adv+" t/dt="+format(t*freq,"04d"))
    fig.savefig('fig1-'+adv+'-'+format(t*freq,"05d")+'.png')
    fig.clf()

  print 'closing the file...'
  nc.close()

  print 'composing animation...'
  cmd = ('convert','*fig1-'+adv+'*.png','fig-'+adv+'.gif')
  subprocess.check_call(cmd)

  print 'deleting frames and the netcdf file...'
  # TODO
