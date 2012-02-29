## @file
#  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
#  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
#  @copyright University of Warsaw
#  @date February 2012
#  @section LICENSE
#    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)

import numpy as np                       # arrays
import subprocess                        # shell calls
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from Scientific.IO.NetCDF import NetCDFFile

nx = 100
nlev = 3
dt = .005
nt = 50
nout = 10

# first: creating a netCDF file with the initial condition
f = NetCDFFile('ini.nc', 'w')

f.createDimension('X', nx)
f.createDimension('Y', 1)
f.createDimension('Z', 1)

v_dp = [None]*nlev
v_qx = [None]*nlev
v_qy = [None]*nlev
for lev in range(nlev) :
  v_dp[lev] = f.createVariable('dp_' + str(lev), 'd', ('X','Y','Z'))
  v_dp[lev][:,0,0] = 1.
  v_qx[lev] = f.createVariable('qx_' + str(lev), 'd', ('X','Y','Z'))
  v_qx[lev][:,0,0] = 0.
  v_qy[lev] = f.createVariable('qy_' + str(lev), 'd', ('X','Y','Z'))
  v_qy[lev][:,0,0] = 0.

v_dp[0][:,0,0] = 1 + .33*pow(np.sin(np.arange(nx) * np.pi / nx),80)
v_dp[1][:,0,0] = 1 - .33*pow(np.sin(np.arange(nx) * np.pi / nx),80)

f.close()

# second: running the model
cmd = (
  '../../icicle',
  '--ini','netcdf',
  '--ini.netcdf.file','ini.nc',
  '--eqs','isentropic',
    '--eqs.isentropic.nlev',str(nlev),
  '--grd.dx','1',
  '--grd.nx',str(nx),
  '--adv','mpdata',
    '--adv.mpdata.fct','1',
    '--adv.mpdata.iord','2',
  '--vel','momeq_extrapol',
  '--nt',str(nt),'--dt',str(dt),'--nout',str(nout),
  '--out','netcdf','--out.netcdf.file','out.nc',
  '--slv','threads','--nsd','1'
)
subprocess.check_call(cmd)

# third: generating a plot
f = NetCDFFile('out.nc', 'r')
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('X [m]')
ax.set_ylabel('time [s]')
ax.set_zlabel('-p [Pa]')

X = f.variables['X']
ny = f.variables['dp_0'].shape[0]
Y = np.arange(0, ny) * nout * dt
X, Y = np.meshgrid(X, Y)

Z = np.zeros((ny, nx)) 
ax.plot_wireframe(X, Y, -Z)

for lev in range(nlev-1,-1,-1) :
  Z += (f.variables['dp_' + str(lev)])[:,:,0,0]
  ax.plot_wireframe(X, Y, -Z)

plt.show()
