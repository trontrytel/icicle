## file
#  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
#  @copyright University of Warsaw
#  @date January 2012
#  @section LICENSE
#    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)

import numpy as np
import subprocess
from scipy.io.netcdf import netcdf_file
from math import sqrt
from sys import exit

file = 'tmp.nc'
sgma = 20
cmd = ( 
  '../../icicle',
  '--ini','gauss',
    '--ini.gauss.A',str(1./sgma/sqrt(2*np.pi)),
    '--ini.gauss.x0','50',
    '--ini.gauss.sx',str(sgma),
  '--grd.dx','1',
  '--grd.nx','100',
  '--adv','upstream',
  '--vel','uniform',
    '--vel.uniform.u','1',
  '--t_max','0',
  '--dt_out','1',
  '--out','netcdf',
    '--out.netcdf.file',file,
    '--out.netcdf.ver','3',
  '--slv','serial'
)
subprocess.check_call(cmd)

nc = netcdf_file(file)
psi = (nc.variables['psi'])[0,:,0,0]
if abs(1.-psi.sum()) > .02 : exit(1)
