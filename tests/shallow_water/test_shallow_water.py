## @file
#  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
#  @copyright University of Warsaw
#  @date January 2012
#  @section LICENSE
#    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)

import numpy as np                       # arrays
import subprocess                        # shell calls
from Scientific.IO.NetCDF import NetCDFFile

nx = 100

# first: creating a netCDF file with the initial condition
f = NetCDFFile('ini.nc', 'w')

f.createDimension('X', nx)
f.createDimension('Y', 1)
f.createDimension('Z', 1)

v_h  = f.createVariable('h', 'd', ('X','Y','Z'))
v_qx = f.createVariable('qx', 'd', ('X','Y','Z'))
v_qy = f.createVariable('qy', 'd', ('X','Y','Z'))

v_h[:,0,0] = 23.
v_qx[:,0,0] = 12.
v_qy[:,0,0] = 13.

f.close()

# second: running the model
cmd = (
  '../../icicle',
  '--ini','netcdf',
  '--ini.netcdf.file','ini.nc',
  '--eqs','shallow_water_2d',
  '--grd.dx','1',
  '--grd.nx',str(nx),
  '--adv','upstream',
  '--vel','uniform','--vel.uniform.u','1', #TODO
  '--t_max','20','--dt_out','1',
  '--out','netcdf','--out.netcdf.file','out.nc',
  '--slv','threads','--nsd','2'
)
subprocess.check_call(cmd)
