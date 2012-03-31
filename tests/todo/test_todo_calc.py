## @file
#  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
#  @copyright University of Warsaw
#  @date March 2012
#  @section LICENSE
#    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)

import numpy as np                       # arrays
import subprocess                        # shell calls
import matplotlib as mpl
import os                                # unlink() etc
import shutil                            # rmtree() etc
from Scientific.IO.NetCDF import NetCDFFile

############################################################################
# simulation parameteters (the 8th WMO Cloud Modelling Workshop: Case 1    #
# W.Grabowski: http://rap.ucar.edu/~gthompsn/workshop2012/case1/case1.pdf  #
############################################################################
nx = 150         # 75     # [1]
ny = 75         # 75     # [1]
dx = 20         # 20     # [m]
dy = 20         # 20     # [m]
z_inv = ny * dy
th_l_0 = 289    # 289    # [K]
q_t_0 = 7.5e-3  # 7.5e-3 # [kg/kg]

############################################################################
# other parameters:
############################################################################
bits = 32
dt = 10.     # [s]
fct = 1  
iord = 2  
nt = 10000
nout = 50 

############################################################################
# physical constants & heper routines for thermodynamics
############################################################################
g = 9.81
cp = 1005.
Rd = 8.314472 / 0.02896
Rdcp = Rd / cp
p0 = 100000 # [Pa] (in definition of theta)

############################################################################
# first: creating a netCDF file with the initial condition                 #
############################################################################
f = NetCDFFile('ini.nc', 'w')

# dimensions
f.createDimension('X', nx)
f.createDimension('Y', ny)
f.createDimension('Z', 1) #TODO: should not be needed

# advected fields
v_rhod_rv = f.createVariable('rhod_rv', 'd', ('X','Y'))

v_rhod_rl = f.createVariable('rhod_rl', 'd', ('X','Y')) # TODO: only for bulk model
v_rhod_rl[:] = 0
v_rhod_rl[:,ny/2:] = .001

v_rhod_th = f.createVariable('rhod_th', 'd', ('X','Y'))

# auxiliary fields
v_rhod = f.createVariable('rhod', 'd', ('X','Y'))

f.close()

############################################################################
# second: running the model                                                #
############################################################################
cmd = (
  '../../icicle',
  '--bits',str(bits),
  '--ini','netcdf',
  '--ini.netcdf.file','ini.nc',
  '--eqs','todo',
  '--grd.dx',str(dx),
  '--grd.nx',str(nx),
  '--grd.dy',str(dy),
  '--grd.ny',str(ny),
  '--adv','mpdata', 
    '--adv.mpdata.fct',str(fct),
    '--adv.mpdata.iord',str(iord),
  '--vel','rasinski',
    '--vel.rasinski.Z_clb',str(z_inv/2),
    '--vel.rasinski.Z_top',str(z_inv),
    '--vel.rasinski.A','1000',
  '--nt',str(nt),'--dt',str(dt),'--nout',str(nout),
  '--out','netcdf','--out.netcdf.file','out.nc',
  '--slv','serial'
)
subprocess.check_call(cmd)
