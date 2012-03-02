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
from Scientific.IO.NetCDF import NetCDFFile

# simulation parameteters (from isopc_yale.::SHAWAT2 by P.K.Smolarkiewicz)
nx = 81 # 160
nz = 10 # 60
dx = 3e3
dz = 600. #300.
dt = 3. # 6.
bv = .012 # [1/s] 
mount_amp = 1.6e3
mount_ro1 = 20e3
uscal = 10 # m/s

nt = 600
nout = 600
#nt = 8000
#nout = 2000

# other parameters 
p_surf = 101300.
th_surf = 300.

# physical constants
g = 9.81
cp = 1005.
Rd = 8.314472 / 0.02896

# initial thermodynamic profile parameters
st = bv**2 / g

# the witch of Agnesi shape
def witch(x, a):
  return 1. / (1 + (x/a)**2)

# potential temperature for constant stability layer
# st = 1/th dth/dz  -->  th = th_0 * exp(st * (z - z0)) 
def theta(z):
  return th_surf * np.exp(st * z)

# pressure profile for hydrostatic constant-stability layer
#         _                                         _  cp/Rd
#        | g p0^Rd/cp  /  1        1  \              |
# p(z) = | ---------  | ------ - ----- | + p_c^Rd/cp |
#        |_  cp st     \ th(z)   th_c /             _|
#                        
def pres(z):
  return (
    g*100000**(Rd/cp)/cp/st*(1/theta(z)-1/th_surf)+p_surf**(Rd/cp)
  )**(cp/Rd)

# first: creating a netCDF file with the initial condition
f = NetCDFFile('ini.nc', 'w')

f.createDimension('X', nx)
f.createDimension('Y', 1) #TODO: should not be needed
f.createDimension('Z', 1) #TODO: should not be needed
f.createDimension('level', nz)

v_dp = [None]*nz
v_qx = [None]*nz
v_qy = [None]*nz
for lev in range(nz) :
  v_dp[lev] = f.createVariable('dp_' + str(lev), 'd', ('X',))
  v_qx[lev] = f.createVariable('qx_' + str(lev), 'd', ('X',))
  v_qy[lev] = f.createVariable('qy_' + str(lev), 'd', ('X',)) # TODO: should not be needed
  v_qy[lev][:] = 0. # TODO: should not be needed

# potential temperatures of the layers (characteristic values)
v_theta = f.createVariable('theta', 'd', ('level',))
for lev in range(nz) :
  v_theta[lev] = theta((lev + .5) * dz)

# topography (temporarily in pressure coordinates)
topo_z = mount_amp * witch(dx * (np.arange(nx+2, dtype='float') - .5 * (nx+1)), mount_ro1)
topo_p = pres(topo_z[1:-1])

# spatial derivatives of the topography
v_dHdx = f.createVariable('dHdx', 'd', ('X',))
v_dHdy = f.createVariable('dHdy', 'd', ('X',)) # should not be needed
v_dHdx[:] = (topo_z[2:] - topo_z[0:-2]) / (2. * dx)
v_dHdy[:] = 0. # TODO: should not be needed

# initilising pressure intervals with topography adjustment
p_up = np.zeros(nx, dtype='float') + pres(nz * dz)
p_dn = np.zeros(nx, dtype='float')
for lev in range(nz-1,-1,-1) :
  for i in range(nx):
    p_dn[i] = min(pres(lev * dz), topo_p[i])
    v_dp[lev][i] = p_dn[i] - p_up[i]
  p_up[:] = p_dn[:]

# initialising momenta
for lev in range(nz-1,-1,-1) :
  for i in range(nx):
    v_qx[lev][i] = uscal * v_dp[lev][i]

# clsing the netCDF tile
f.close()

# second: running the model
cmd = (
  '../../icicle',
  '--ini','netcdf',
  '--ini.netcdf.file','ini.nc',
  '--eqs','isentropic',
    '--eqs.isentropic.nlev',str(nz),
    '--eqs.isentropic.p_max',str(pres(nz * dz)),
  '--grd.dx',str(dx),
  '--grd.nx',str(nx),
  '--adv','mpdata',
    '--adv.mpdata.fct','1',
    '--adv.mpdata.iord','2',
  '--vel','momeq_extrapol',
  '--nt',str(nt),'--dt',str(dt),'--nout',str(nout),
  '--out','netcdf','--out.netcdf.file','out.nc',
  '--slv','serial','--nsd','1'
)
subprocess.check_call(cmd)

# third: generating a plot
f = NetCDFFile('out.nc', 'r')
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlabel('X [m]')
ax.set_ylabel('-ln(p/p_surf) [1]')

for t in range(nt/nout+1):
  X = f.variables['X']

  # top isentrope
  P = np.zeros(nx) + pres(nz * dz)
  ax.plot(X, -np.log(P/p_surf))

  # other isentropes
  for lev in range(nz-1,-1,-1) :
    P += (f.variables['dp_' + str(lev)])[t,:,0,0]
    ax.plot(X, -np.log(P/p_surf))

ax.plot(X, -np.log(topo_p/p_surf))
plt.savefig('fig1.svg')
