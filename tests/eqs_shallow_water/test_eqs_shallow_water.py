## @file
#  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
#  @copyright University of Warsaw
#  @date January 2012 - February 2012
#  @section LICENSE
#    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)

import numpy as np                       # arrays
import subprocess                        # shell calls
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from Scientific.IO.NetCDF import NetCDFFile

nx = 20
dt = .075
nt = 100
nout = 5

# first: creating a netCDF file with the initial condition
f = NetCDFFile('ini.nc', 'w')

f.createDimension('X', nx)
f.createDimension('Y', 1)
f.createDimension('Z', 1)

v_h  = f.createVariable('h', 'd', ('X','Y','Z'))
v_qx = f.createVariable('qx', 'd', ('X','Y','Z'))
v_qy = f.createVariable('qy', 'd', ('X','Y','Z'))

v_h[:,0,0] = 20.
v_h[nx/3:nx/2,0,0] = 20.1
v_qx[:,0,0] = 0.
v_qy[:,0,0] = 0.

f.close()

# second: running the model
cmd = (
  '../../icicle',
  '--ini','netcdf',
  '--ini.netcdf.file','ini.nc',
  '--eqs','shallow_water',
  '--grd.dx','1',
  '--grd.nx',str(nx),
  '--adv','mpdata',
    '--adv.mpdata.fct','0',
    '--adv.mpdata.iord','2',
  '--vel','momeq_extrapol',
  '--nt',str(nt),'--dt',str(dt),'--nout',str(nout),
  '--out','netcdf','--out.netcdf.file','out.nc',
  '--slv','serial'
)
subprocess.check_call(cmd)

# third: generating a plot
f = NetCDFFile('out.nc', 'r')
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
X = f.variables['X']
Y = np.arange(0, f.variables['h'].shape[0]) #* getattr(f, dt_out)
Z = (f.variables['h'])[:,:,0,0]
print X.shape, Y.shape, Z.shape
X, Y = np.meshgrid(X, Y)
ax.plot_wireframe(X, Y, Z)
plt.show()
