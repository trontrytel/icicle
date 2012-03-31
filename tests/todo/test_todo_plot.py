## @file
#  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
#  @copyright University of Warsaw
#  @date March 2012
#  @section LICENSE
#    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)

import numpy as np                       # arrays
import subprocess                        # shell calls
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from Scientific.IO.NetCDF import NetCDFFile
import os                                # unlink() etc
import shutil                            # rmtree() etc

############################################################################
# generating a plot
############################################################################
f = NetCDFFile('out.nc', 'r')
# nt = f.sync() # does not work
# nt = f.dimensions['time'] # does not work
nt = f.variables['rhod_th'].shape[0]
dt_out = f.variables['dt_out'].getValue()
X = f.variables['X']
Y = f.variables['Y']

X,Y = np.meshgrid(X,Y) 
X = np.transpose(X)
Y = np.transpose(Y)

os.mkdir('tmp')

fig = plt.figure()
for t in range(nt):
  plt.clf()
  ax = fig.add_subplot(111)
  ax.set_xlabel('X [m]')
  ax.set_ylabel('Y [m]')
  plt.suptitle('t = ' + format(int(t * dt_out), '05d') + ' s')
  
  V = f.variables['rhod_rl'][t,:,:,0]
  print V.shape
  contour = plt.contourf(X,Y, 1e3 * V)
  cb = plt.colorbar(contour)
  cb.set_label("liquid water mass density [g/kg]")
  plt.savefig('tmp/frame_'+format(t,"05d")+'.png', dpi=125)

f.close
cmd=('convert','tmp/frame_*.png','todo.gif')
subprocess.check_call(cmd)
shutil.rmtree('tmp')
