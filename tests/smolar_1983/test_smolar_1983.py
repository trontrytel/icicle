import numpy as np                       # arrays
from scipy.io.netcdf import netcdf_file  # netcdf IO
import subprocess                        # shell calls
import os                                # unlink()
import matplotlib                        # plotting
matplotlib.use('Agg')                    # non-interactive backend
import matplotlib.pyplot as plt          # ?
from mpl_toolkits.mplot3d import Axes3D  # 3D plot
from matplotlib import cm                # colormap 

print 'running the model...'
if os.path.exists('test_smolar_1983.nc') : os.unlink('test_smolar_1983.nc')
cmd = (
  '../../icicles',
  '--dt','.1', 
  '--grd.dx','1',
  '--grd.dy','1',
  '--grd.dz','1',
  '--nt','628', # should be 6*628 !
  '--nx','100',
  '--nz','100',
  '--vel','test',
  '--vel.test.omega','.1',
  '--adv','leapfrog',
  '--ini','cone', 
  '--ini.cone.h','3.87', 
  '--ini.cone.x0','75', 
  '--ini.cone.z0','50', 
  '--ini.cone.r','15',
  '--slv','threads','--nsd','2',
  '--out','netcdf', 
  '--out.netcdf.ver','3', 
  '--out.netcdf.freq','10', 
  '--out.netcdf.file','test_smolar_1983.nc'
)
subprocess.check_call(cmd)

print 'opening the output file...' # TODO: move it into the loop!
nc = netcdf_file("test_smolar_1983.nc")
psi = nc.variables['psi']

print 'setting-up the plot...'
x = np.arange(.5, 100.5, 1)
y = np.arange(.5, 100.5, 1)
x, y = np.meshgrid(x, y)

print 'saving animation frames...'
fig = plt.figure()
for t in range(psi.shape[0]) :
  axs = fig.gca(projection='3d')
  srf = axs.plot_surface(x, y, psi[t,:,0,:], rstride=1, cstride=1, cmap=cm.jet, linewidth=0, antialiased=True, vmin=0, vmax=4)
  axs.set_zlim3d(0., 4.)
  fig.colorbar(srf) #shrink=0.5, aspect=5)
  fig.savefig('fig1-'+format(t,"05d")+'.png')
  fig.clf()

print 'closing the file...'
nc.close()

print 'composing animation...'
cmd = ('convert','*fig1-*.png','fig.gif')
subprocess.check_call(cmd)

print 'deleting frames and the netcdf file...'
# TODO
