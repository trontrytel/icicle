## @file
#  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
#  @copyright University of Warsaw
#  @date Januar 2012
#  @section LICENSE
#    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)

#  1d isolines test from pks & wwg 1990 

import numpy as np                       # arrays
from scipy.io.netcdf import netcdf_file  # netcdf IO
import subprocess                        # shell calls
import os                                # unlink()
import sys                               # exit
import math                              # sqrt
import matplotlib.pyplot as plt          # plots
from matplotlib.mlab import griddata     # griddata

#dx, cour, err = np.loadtxt('new_iord2_fct0_toa1_b64.txt', unpack=True)
dx, cour, err = np.loadtxt('isolines.txt', unpack=True)

theta=np.zeros(dx.shape[0])
r=np.zeros(dx.shape[0])
x=np.zeros(dx.shape[0])
y=np.zeros(dx.shape[0])

norm=max(dx)
theta=cour*math.pi/2.
for i in range(err.shape[0]) :
  err[i]=math.log(err[i],2)
for i in range(dx.shape[0]) :
  r[i]=math.log(dx[i]/norm,2)+8
for i in range(theta.shape[0]) :
  x[i]=r[i]*math.cos(theta[i])
  y[i]=r[i]*math.sin(theta[i])

ngrid = 800*2

xi = np.linspace(0, 8, ngrid)
yi = np.linspace(0, 8, ngrid)
zi = griddata(x,y,err,xi,yi,interp='linear')

fig = plt.figure()
plt.contour(xi,yi,zi,15,linewidths=0.5,colors='k')
plt.contourf(xi,yi,zi,15,cmap=plt.cm.jet)
plt.colorbar() # draw colorbar

plt.xlim(0,8)
plt.ylim(0,8)
plt.title('log2(err)')
plt.xlabel('r; C=0')
plt.ylabel('r; C=1')

#plt.show()
fig.savefig('new_iord2_fct1_toa0_b64.pdf')
fig.clf
