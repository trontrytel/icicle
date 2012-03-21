## @file
#  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
#  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
#  @copyright University of Warsaw
#  @date February 2012
#  @section LICENSE
#    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)

import numpy as np                       # arrays
import subprocess                        # shell calls
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import MaxNLocator
from Scientific.IO.NetCDF import NetCDFFile
import os                                # unlink() etc
import shutil                            # rmtree() etc
import scitools
import math

############################################################################
# simulation parameteters (from isopc_yale.::SHAWAT2 by P.K.Smolarkiewicz) #
############################################################################
bits = 32
nx = 300     # 161  # [1]
nz = 60      # 61   # [1]
dx = 3e3     # 3e3  # [m]
dz = 300.    # 300. # [m]
dt = 6.      # 6.   # [s]
bv = .012    # .012 # [1/s] (Brunt-Vaisala frequency)
mount_amp = 1.6e3 # 1.6e3
mount_ro1 = 20e3  # 20e3
abslev = int(3. * nz / 4.) # gravity-wave absorber starts at .75 of the domain height
absamp = 1 # 1 ?
uscal = 10.  # 10. # [m/s]
p_surf = 101300.
th_surf = 300.
fct = 1 # 1
iord = 2 # 2
nt = 1000 # 8000
nout = 25 # 2000

############################################################################
# physical constants & heper routines for thermodynamics
############################################################################
g = 9.81
cp = 1005.
Rd = 8.314472 / 0.02896
Rdcp = Rd / cp
p0 = 100000 # [Pa] (in definition of theta)
st0 = bv**2 / g

# the witch of Agnesi shape
def witch(x, a):
  return 1. / (1 + (x/a)**2)

# potential temperature for constant stability layer
# st = 1/th dth/dz  -->  th = th_0 * exp(st * (z - z0)) 
def theta(z, st, th_c):
  return th_c * np.exp(st * z)

# pressure profile for hydrostatic constant-stability layer
#         _                                         _  cp/Rd
#        | g p0^Rd/cp  /  1        1  \              |
# p(z) = | ---------  | ------ - ----- | + p_c^Rd/cp |
#        |_  cp st     \ th(z)   th_c /             _|
#                        
def pres(z, st, th_c, p_c):
  return (
    g*p0**(Rdcp)/cp/st*(1/theta(z,st,th_c)-1/th_c)+p_c**(Rdcp)
  )**(1./Rdcp)

############################################################################
# first: creating a netCDF file with the initial condition                 #
############################################################################
f = NetCDFFile('ini.nc', 'w')

f.createDimension('X', nx)
f.createDimension('Y', 1) #TODO: should not be needed
f.createDimension('Z', 1) #TODO: should not be needed
f.createDimension('dlevel', nz - 1)

v_dp = [None]*nz
v_qx = [None]*nz
v_qy = [None]*nz
for lev in range(nz) :
  v_dp[lev] = f.createVariable('dp_' + str(lev), 'd', ('X',))
  v_qx[lev] = f.createVariable('qx_' + str(lev), 'd', ('X',))

# potential temperatures of the layers (characteristic values)
v_dtheta = f.createVariable('dtheta', 'd', ('dlevel',))
for lev in range(nz-1) :
  v_dtheta[lev] = theta((lev+1.5) * dz, st0, th_surf) - theta((lev+.5) * dz, st0, th_surf) 

theta_frst = theta(.5 * dz, st0, th_surf)

# topography and its spatial derivatives of the topography
topo_z = mount_amp * witch(dx * (np.arange(nx+2, dtype='float') - .5 * (nx+1)), mount_ro1)
v_dHdx = f.createVariable('dHdx', 'd', ('X',))
v_dHdx[:] = (topo_z[2:] - topo_z[0:-2]) / (2. * dx)
topo_z = topo_z[1:-1]

# initilising pressure intervals with topography adjustment
topo_p = pres(topo_z, st0, th_surf, p_surf)
p_top = pres(nz * dz, st0, th_surf, p_surf)
p_up = np.zeros(nx, dtype='float') + p_top
p_dn = np.zeros(nx, dtype='float')
for lev in range(nz-1,-1,-1) :
  for i in range(nx):
    p_dn[i] = min(pres(lev * dz, st0, th_surf, p_surf), topo_p[i])
    v_dp[lev][i] = p_dn[i] - p_up[i]
  p_up[:] = p_dn[:]

# saving density
rho = np.zeros(nz, dtype='float')
for lev in range(nz) :
  rho[lev] = v_dp[lev][0] / g / dz

# initialising momenta
for lev in range(nz-1,-1,-1) :
  for i in range(nx):
    v_qx[lev][i] = uscal * v_dp[lev][i]

# clsing the netCDF tile
f.close()

############################################################################
# second: running the model                                                #
############################################################################
cmd = (
  '../../icicle',
  '--bits',str(bits),
  '--ini','netcdf',
  '--ini.netcdf.file','ini.nc',
  '--eqs','isentropic',
    '--eqs.isentropic.nlev',str(nz),
    '--eqs.isentropic.p_top',str(p_top),
    '--eqs.isentropic.theta_frst',str(theta_frst),
    '--eqs.isentropic.abslev',str(abslev),
    '--eqs.isentropic.absamp',str(absamp),
  '--grd.dx',str(dx),
  '--grd.nx',str(nx),
  '--adv','mpdata', 
    '--adv.mpdata.fct',str(fct),
    '--adv.mpdata.iord',str(iord),
  '--vel','momeq_extrapol',
  '--nt',str(nt),'--dt',str(dt),'--nout',str(nout),
  '--out','netcdf','--out.netcdf.file','out.nc',
  '--slv','serial'
)
subprocess.check_call(cmd)

############################################################################
# third: generating a plot
############################################################################
f = NetCDFFile('out.nc', 'r')
os.mkdir('tmp')

qui_X = np.zeros((nx,nz))
qui_Z = np.zeros((nx,nz))
qui_U = np.zeros((nx,nz))
qui_V = np.zeros((nx,nz))
qui_C = np.zeros((nx,nz))
for lev in range(nz): qui_X[:,lev] = f.variables['X'] 
qui_X *= 1e-3 # m -> km

fig = plt.figure()
for t in range(nt/nout+1):
  plt.clf()
  ax = fig.add_subplot(111)
  ax.set_xlabel('X [km]')
  ax.xaxis.set_major_locator(MaxNLocator(4))
  plt.ylim([0,.9*nz*dz])
  dom_half_km = .5 * nx * dx * 1e-3
  mount_ro1_km = mount_ro1 * 1e-3
  plt.xlim([dom_half_km - 3 * mount_ro1_km, dom_half_km + 3 * mount_ro1_km])
  ax.set_ylabel('h [m]')
  ax.fill(qui_X[:,0], topo_z, color='#BBBBBB', linewidth=0)
  plt.suptitle('t = ' + format(int(t * dt), '05d') + ' s')
  
  # from top to bottom - calculating pressures
  p_dn = np.zeros(nx) + p_top
  for lev in range(nz-1,-1,-1) :
    p_dn[:] += (f.variables['dp_' + str(lev)])[t,:,0,0]

  # from bottom to top - calculating heights
  z_dn = topo_z
  z_up = np.zeros(nx)
  for lev in range(nz) :
    # calculating isentrope height
    dp = (f.variables['dp_' + str(lev)])[t,:,0,0]
    qx = (f.variables['qx_' + str(lev)])[t,:,0,0]
    z_up = z_dn + dp / (rho[lev] * g)
    # registering values for the plot
    qui_Z[:,lev] = z_dn+(z_up-z_dn)/2. 
    for i in range(nx) :
      # only for cyclic conditions!
      l = (i-1)%nx
      r = (i+1)%nx
      qui_U[i,lev] = 2*dx / math.pow( math.pow(z_dn[r]-z_dn[l],2) + math.pow(2*dx,2) , .5) * qx[i] / dp[i]
      qui_V[i,lev] = (z_dn[r]-z_dn[l]) / math.pow( math.pow(z_dn[r]-z_dn[l],2) + math.pow(2*dx,2) , .5) * math.fabs(qx[i]) / dp[i]
      qui_C[i,lev] = math.pow(math.pow(qui_U[i,lev],2) + math.pow(qui_V[i,lev],2),.5) 
    # storing values for the next iteration
    z_dn = z_up
    p_dn -= dp

  vectors = ax.quiver(qui_X, qui_Z, qui_U/dx, qui_V/dz, qui_C, 
    units='dots', pivot='mid', width=1, headwidth=4.5, headlength=6, headaxislength=6
  )
  cb = plt.colorbar(vectors)
  cb.set_label("air velocity [m/s]")
  vectors.set_clim([0,25])
  plt.savefig('tmp/frame_'+format(t,"05d")+'.png', dpi=125)

f.close
cmd=('convert','tmp/frame_*.png','halny_vel.gif')
subprocess.check_call(cmd)
shutil.rmtree('tmp')
