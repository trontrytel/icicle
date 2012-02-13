## @file
#  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
#  @copyright University of Warsaw
#  @date January 2012 - February 2012
#  @section LICENSE
#    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)

import numpy as np                       # arrays
import subprocess                        # shell calls
from Scientific.IO.NetCDF import NetCDFFile

nx = 5

# first: creating a netCDF file with the initial condition
f = NetCDFFile('ini.nc', 'w')

f.createDimension('X', nx)
f.createDimension('Y', 1)
f.createDimension('Z', 1)

v_h  = f.createVariable('h', 'd', ('X','Y','Z'))
v_qx = f.createVariable('qx', 'd', ('X','Y','Z'))
v_qy = f.createVariable('qy', 'd', ('X','Y','Z'))

v_h[:,0,0] = 23.
v_h[0,0,0] = 24.
v_qx[:,0,0] = 1
v_qy[:,0,0] = 13.

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
    '--adv.mpdata.fct','0',
    '--adv.mpdata.iord','3',
  '--vel','momeq_extrapol',
  '--nt','20','--dt','1','--nout','1',
  '--out','netcdf','--out.netcdf.file','out.nc',
  '--slv','serial'
)
subprocess.check_call(cmd)
