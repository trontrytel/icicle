## @file
#  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
#  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
#  @copyright University of Warsaw
#  @date March 2012
#  @section LICENSE
#    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)

import numpy as np                       # arrays
import subprocess                        # shell calls
import matplotlib as mpl 
mpl.use('Agg') # non-X terminal
import matplotlib.pyplot as plt 
from Scientific.IO.NetCDF import NetCDFFile
import os                                # unlink() etc
import shutil                            # rmtree() etc

nx = 1000
dx = 1.
dt = 1.
nt = 1600
nout = 10
omega =2*np.pi / dt / 400
u = 0.5

# first: creating a netCDF file with the initial condition
f = NetCDFFile('ini.nc', 'w')
f.createDimension('X', nx)
f.createDimension('Y', 1) # TODO
f.createDimension('Z', 1) # TODO
v_psi = f.createVariable('psi', 'd', ('X',))
v_phi = f.createVariable('phi', 'd', ('X',))
v_psi[:] = pow(np.sin( (np.arange(nx) * np.pi / nx) ), 300)
v_psi[:] = np.roll(v_psi[:] , -nx/2+80)
v_phi[:] = 0
f.close()

# second: running the model
cmd = (
  '../../icicle',
  '--ini','netcdf',
  '--ini.netcdf.file','ini.nc',
  '--eqs','harmonic_oscillator',
    '--eqs.harmonic_oscillator.omega',str(omega),
  '--grd.dx',str(dx),
  '--grd.nx',str(nx),
  '--adv','mpdata',
    '--adv.mpdata.fct','1',
    '--adv.mpdata.iord','2',
  '--vel','uniform',
    '--vel.uniform.u',str(u),
  '--nt',str(nt),'--dt',str(dt),'--nout',str(nout),
  '--out','netcdf','--out.netcdf.file','out.nc',
  '--slv','threads','--nsd','2'
)
subprocess.check_call(cmd)

# third: generating a plot
f = NetCDFFile('out.nc', 'r')
os.mkdir('tmp')
for t in range(nt/nout):
  fig = plt.figure()
  plt.ylim((-1, 1))
  ax = fig.add_subplot(111)
  ax.set_xlabel('X')
  X = f.variables['X']
  psi = (f.variables['psi'])[t,:,0,0]
  phi = (f.variables['phi'])[t,:,0,0]
  ax.plot(X, psi, color='#FF0000')
  ax.plot(X, phi, color='#0000FF')
  plt.savefig('tmp/frame_'+format(t,"05d")+'.png')

f.close
cmd=('convert','tmp/frame_*.png','harmonic_osc.gif')
subprocess.check_call(cmd)
shutil.rmtree('tmp')


