## @file
#  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
#  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
#  @copyright University of Warsaw
#  @date Januar 2012
#  @section LICENSE
#    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)

#  1d isolines test from pks & wwg 1990 

import numpy as np                       # arrays
from scipy.io.netcdf import netcdf_file  # netcdf IO
import subprocess                        # shell calls
import os                                # unlink()
import sys                               # exit
import math                              # sqrt
import matplotlib.pyplot as plt          # plots
from matplotlib.ticker import AutoMinorLocator

file = 'tmp.nc'

iord = 2 
fct = 0
toa = 0
bits = 32
print "running 1D MPDATA simultation with iord=", str(iord), " fct=", str(fct), " toa=", str(toa)

f = open('isolines.txt','w')
# simulation parameters
t_max = 1. # "arbitrarily"
dx_max = 1.
x_max = 20.*dx_max
velocity = dx_max / t_max # "solution advects over the one grid increment for r=8"
sgma = 1.5*dx_max
x0 = .5*x_max

for dx in (dx_max, dx_max/2, dx_max/4, dx_max/8, dx_max/16, dx_max/32, dx_max/64, dx_max/128) : 
  for cour in (.05, .1, .15, .2, .25, .3, .35, .4, .45, .5, .55, .6, .65, .7, .75, .8, .85, .9, .95) :
    dt = cour*dx/velocity
    nx = int(round(x_max/dx))
    nt = int(round(t_max/dt))
    print "dt", dt
    print "dx", dx
    print "courant", cour
    print "nt", nt
    print "nt*dt", nt*dt
    cmd = (
      '../../icicle',
      '--bits',str(bits),
      '--dt',str(dt), 
      '--grd.dx',str(dx),
      '--nt',str(nt),
      '--nout',str(nt),
      '--grd.nx',str(nx),
      '--vel','uniform',
      '--vel.uniform.u',str(velocity),
      '--adv','mpdata',
        '--adv.mpdata.fct',str(fct),
        '--adv.mpdata.third_order',str(toa),
        '--adv.mpdata.iord',str(iord),
      '--ini','gauss',
        '--ini.gauss.A',str(1./sgma/math.sqrt(2*np.pi)),
        '--ini.gauss.x0',str(x0),
        '--ini.gauss.sx',str(sgma),
      '--slv', 'serial',
      '--out','netcdf', 
      '--out.netcdf.ver','3', 
      '--out.netcdf.file',file
      )
    if os.path.exists(file) : os.unlink(file)
    subprocess.check_call(cmd)

    #read from netcdf file  
    nc = netcdf_file(file)
    location = (nc.variables['X'])[:]
    x_k=x0+velocity*dt*nt
    exact=np.zeros(location.shape[0])
    for i in range(location.shape[0]) :
      exact[i]=1./sgma/math.sqrt(2*np.pi)*math.exp(-.5*math.pow(location[i]-x_k,2)*math.pow(sgma,-2))
    psint = (nc.variables['psi'])[1,:,0,0]
    time = float(nc.user_time_seconds)
    err=math.sqrt(((psint - exact)*(psint - exact)).sum()/nx)/(nt * dt)

    #cleanup
    print 'closing and deleting the file...'
    nc.close()
    if os.path.exists(file) : os.unlink(file)

    f.write(str(dx) + ' ' + str(cour) + ' ' + str(err) + '\n')
f.close()
