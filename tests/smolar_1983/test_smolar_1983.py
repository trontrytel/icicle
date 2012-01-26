## @file
#  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
#  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
#  @copyright University of Warsaw
#  @date December 2011
#  @section LICENSE
#    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)

import numpy as np                       # arrays
from scipy.io.netcdf import netcdf_file  # netcdf IO
import subprocess                        # shell calls
import os                                # unlink() etc
import shutil                            # rmtree() etc
import matplotlib                        # plotting
matplotlib.use('Agg')                    # non-interactive backend
import matplotlib.pyplot as plt          # ?
from mpl_toolkits.mplot3d import Axes3D  # 3D plot
from matplotlib import cm                # colormap 
from matplotlib.colors import LinearSegmentedColormap 

for adv in ('mpdata',):#'leapfrog','upstream'):
  print 'running the model...'
  os.mkdir('tmp')
  file = 'tmp/test_smolar_1983'+adv+'.nc'
  if os.path.exists(file) : os.unlink(file)
  freq = 25
  dt = .1
  cmd = (
    '../../icicle',
    '--dt',str(dt), 
    '--grd.dx','1',
    '--grd.dy','1',
    '--grd.dz','1',
    '--t_max','188.4', # 6*628 in the original paper !
    '--grd.nx','100',
    '--grd.nz','100',
    '--vel','test',
    '--vel.test.omega','.1',
    '--adv',adv,
      '--adv.mpdata.fct','1',
      '--adv.mpdata.iord','3',
      '--adv.mpdata.third_order','1',
    '--ini','cone', 
    '--ini.cone.h','3.87', 
    '--ini.cone.x0','75', 
    '--ini.cone.z0','50', 
    '--ini.cone.r','15',
    '--ini.cone.h0','1',
    '--slv','threads','--nsd','2',
    '--out','netcdf', 
    '--out.netcdf.ver','3', 
    '--dt_out',str(freq * dt),
    '--out.netcdf.file',file
  )
  subprocess.check_call(cmd)

  print 'opening the output file...' 
  nc = netcdf_file(file)
  psi = nc.variables['psi']

  print 'setting-up the plot...'
  x = np.arange(.5, 100.5, 1)
  y = np.arange(.5, 100.5, 1)
  x, y = np.meshgrid(x, y)

  print 'saving animation frames...'
  vmin = 0
  vmax = 5.
  mycmdata2 = {
    'red' :  ((0., 1., 1.), (0.1999, .5, 0.), (1., 1., 1.)),
    'green': ((0., 0., 0.), (0.1999, 0., .1), (1., 1., 1.)), 
    'blue' : ((0., 0., 0.), (0.1999, 0., 1.), (.6, 0., 0.), (1., 0., 0.))
  }
  mycm = LinearSegmentedColormap('mycm', mycmdata2) 
  fig = plt.figure()
  for t in range(psi.shape[0]) :
    axs = fig.gca(projection='3d')
    srf = axs.plot_surface(x, y, psi[t,:,0,:], rstride=1, cstride=1, cmap=mycm, linewidth=0, antialiased=True, vmin=vmin, vmax=vmax)
    axs.set_zlim3d(vmin, vmax)
    fig.colorbar(srf)
    fig.suptitle(adv+" t/dt="+format(t*freq,"04d"))
    fig.savefig('tmp/'+adv+'-'+format(t*freq,"05d")+'.png')
    fig.clf()

  print 'closing the file...'
  nc.close()

  print 'composing animation...'
  cmd = ('convert','tmp/'+adv+'*.png','fig-'+adv+'-2d.gif')
  subprocess.check_call(cmd)

  print 'deleting temporary file (animation frames and the netcdf file)...'
  shutil.rmtree('tmp')
