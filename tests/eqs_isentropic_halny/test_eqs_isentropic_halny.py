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
from Scientific.IO.NetCDF import NetCDFFile

############################################################################
# simulation parameteters (from isopc_yale.::SHAWAT2 by P.K.Smolarkiewicz) #
############################################################################
nx = 160     # 161  # [1]
nz = 60      # 61   # [1]
dx = 3e3     # 3e3  # [m]
dz = 300.    # 300. # [m]
dt = 6.      # 6.   # [s]
bv = .012    # .012 # [1/s] (Brunt-Vaisala frequency)
mount_amp = 1.6e3 # 1.6e3
mount_ro1 = 20e3  # 20e3
abslev = nz +1 # 3 * (nz / 4) # gravity-wave absorber starts at .75 of the domain height
absamp = 1
uscal = 10.  # 10. # [m/s]
p_surf = 101300.
th_surf = 300.
fct = 1
iord = 2
nt = 50 # 8000
nout = nt # 2000

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

# height of the top of a layer
def alt(th_up, th_dn, st, z_dn):
  return z_dn + 1./st * np.log(th_up / th_dn)

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

# stability within a constant-stability layer defined by 
# thetas and pressures at the bottom and top
#       g       p0^Rd/cp        / 1     1   \
# st = --- ------------------- | --- - ----  |
#      cp  p^Rd/cp - p_c^Rd/cp  \ th   th_c /
def sta(p_dn, p_up, th_dn, th_up):
  return g / cp * (1/th_up - 1/th_dn) * p0**(Rdcp) / (p_up**Rdcp - p_dn**Rdcp)

############################################################################
# first: creating a netCDF file with the initial condition                 #
############################################################################
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
  #v_qy[lev] = f.createVariable('qy_' + str(lev), 'd', ('X',)) # TODO: should not be needed
  #v_qy[lev][:] = 0. # TODO: should not be needed

# potential temperatures of the layers (characteristic values)
v_dtheta = f.createVariable('dtheta', 'd', ('level',))
dtheta = np.zeros(nz, dtype='float')
for lev in range(nz) :
  dtheta[lev] = theta((lev+1) * dz, st0, th_surf) - theta((lev) * dz, st0, th_surf) 
  v_dtheta[lev] = dtheta[lev]

# topography and its spatial derivatives of the topography
topo_z = mount_amp * witch(dx * (np.arange(nx+2, dtype='float') - .5 * (nx+1)), mount_ro1)
v_dHdx = f.createVariable('dHdx', 'd', ('X',))
v_dHdx[:] = (topo_z[2:] - topo_z[0:-2]) / (2. * dx)
v_dHdy = f.createVariable('dHdy', 'd', ('X',)) # should not be needed
v_dHdy[:] = 0. # TODO: should not be needed
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
  '--ini','netcdf',
  '--ini.netcdf.file','ini.nc',
  '--eqs','isentropic',
    '--eqs.isentropic.nlev',str(nz),
    '--eqs.isentropic.p_max',str(p_top),
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
  '--slv','serial','--nsd','1'
)
subprocess.check_call(cmd)

# third: generating a plot
f = NetCDFFile('out.nc', 'r')
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlabel('X [m]')
ax.set_ylabel('h [m]')

X = f.variables['X']
ax.fill(X, topo_z, color='#BBBBBB', linewidth=0)
#plt.ylim([0,nz*dz])
plt.xlim([0,nx*dx])

# TODO: animation!
for t in range(nt/nout+1):
  # from top to bottom - calculating pressures
  p_dn = np.zeros(nx) + p_top
  for lev in range(nz-1,-1,-1) :
    p_dn[:] += (f.variables['dp_' + str(lev)])[t,:,0,0]

  # from bottom to top - calculating heights
  th_dn = th_surf
  z_dn = topo_z
  z_up = np.zeros(nx)
  for lev in range(nz) :
    # calculating isentrope height
    dp = (f.variables['dp_' + str(lev)])[t,:,0,0]
    th_up = th_dn + dtheta[lev]
    st = sta(p_dn, p_dn - dp, th_dn, th_up)
    z_up = alt(th_up, th_dn, st, z_dn)
    # plotting
    ax.plot(X, z_dn, color=cm.gist_heat((1.*lev)/nz,1))
    ax.plot(X, z_up, color=cm.gist_heat((1.*lev)/nz,1))
    # storing values for the next iteration
    z_dn = z_up
    th_dn = th_up
    p_dn -= dp

plt.savefig('fig1.svg')
