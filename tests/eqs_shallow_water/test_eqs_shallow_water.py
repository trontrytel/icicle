## @file
#  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
#  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
#  @copyright University of Warsaw
#  @date January 2012 - February 2012
#  @section LICENSE
#    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)

import numpy as np                       # arrays
import subprocess                        # shell calls
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from Scientific.IO.NetCDF import NetCDFFile

nx = 100
dt = .05
nt = 250
nout = 10

# first: creating a netCDF file with the initial condition
f = NetCDFFile('ini.nc', 'w')

f.createDimension('X', nx)
f.createDimension('Y', 1) # TODO: should not be needed
f.createDimension('Z', 1) # TODO: should not be needed

v_h  = f.createVariable('h', 'd', ('X',))
v_qx = f.createVariable('qx', 'd', ('X',))
v_qy = f.createVariable('qy', 'd', ('X',)) # TODO: should not be needed
v_dHdx = f.createVariable('dHdx', 'd', ('X',))
v_dHdy = f.createVariable('dHdy', 'd', ('X',)) # TODO: should not be needed

v_h[:] = 1. - .1*pow(np.sin(np.arange(nx) * np.pi / nx),80)
v_qx[:] = 0.
v_qy[:] = 0.
v_dHdx[:] = 0.
v_dHdy[:] = 0.

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
    '--adv.mpdata.fct','1',
    '--adv.mpdata.iord','2',
  '--vel','momeq_extrapol',
  '--nt',str(nt),'--dt',str(dt),'--nout',str(nout),
  '--out','netcdf','--out.netcdf.file','out.nc',
  '--slv','threads','--nsd','2'
)
subprocess.check_call(cmd)

# third: generating a plot
f = NetCDFFile('out.nc', 'r')
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('X [m]')
ax.set_ylabel('time [s]')
ax.set_zlabel('h [m]')
X = f.variables['X']
Y = np.arange(0, f.variables['h'].shape[0]) * nout * dt
Z = (f.variables['h'])[:,:,0,0]
X, Y = np.meshgrid(X, Y)
ax.plot_wireframe(X, Y, Z)
plt.savefig('fig1.svg')
