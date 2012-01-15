## @file
#  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
#  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
#  @copyright University of Warsaw
#  @date Januar 2012
#  @section LICENSE
#    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)

import numpy as np                       # arrays
from scipy.io.netcdf import netcdf_file  # netcdf IO
import subprocess                        # shell calls
import os                                # unlink()
import sys                               # exit
import math                              # sqrt
from texttable import Texttable          # ascii-art tables
import matplotlib.pyplot as plt          # plots

table = Texttable()
table.set_deco(Texttable.HEADER + Texttable.VLINES)
table.add_rows([
  ["iord","fct","toa","max-1.1","min-0.1","L1_err","L2_err","cpu [s]"]
])
table.set_cols_width(
  [ 4,     3,    3,    7,        7,        7,       7,       8]
)

file = 'tmp.nc'
time_ref = 1.
fig = 0
for iord in (1,2,3,4) : #,40) : 
  for fct in (0,1):
    if iord == 1 and fct == 1 : continue # FCT works from iord > 1
    for toa in (0,1):
      if iord == 1 and toa == 1 : continue # third_order opt. works from iord > 1
      if iord == 40 and (fct == 1 or toa == 1) : continue # just "BASIC MPDATA"
      for bits in (32,) : # ,64,128) : 
        fig += 1
        print "running 1D MPDATA simultation with iord=", str(iord), " fct=", str(fct), " toa=", str(toa)
    
        # simulation parameters
        nt = 40
        courant = .5
        cmd = (
          '../../icicle',
          '--bits',str(bits),
          '--dt','1', 
          '--grd.dx','1',
          '--nt',str(nt),
          '--nx',str(60),
          '--vel','uniform',
          '--vel.uniform.u',str(courant),
          '--adv','mpdata',
            '--adv.mpdata.fct',str(fct),
            '--adv.mpdata.third_order',str(toa),
            '--adv.mpdata.iord',str(iord),
          '--ini','boxcar', 
          '--ini.boxcar.ax','0', 
          '--ini.boxcar.bx','20', 
          '--ini.boxcar.A','1', 
          '--ini.boxcar.A0','.1', 
          '--slv','serial',
          '--out','netcdf', 
          '--out.netcdf.ver','3', 
          '--out.netcdf.freq',str(nt), 
          '--out.netcdf.file',file
        )
        if os.path.exists(file) : os.unlink(file)
        subprocess.check_call(cmd)
  
        nc = netcdf_file(file)
        exact = (nc.variables['psi'])[0,:,0,0]
        exact = np.roll(exact, int(courant * nt), 0)
        psint = (nc.variables['psi'])[1,:,0,0]
        time = float(nc.user_time_seconds)
        #if bits == 32 and iord == 1 : time_ref = time
  
        plt.subplot(4,4,fig)
        plt.annotate("iord=" + str(iord) + "\nfct=" + str(fct) + "\ntoa=" + str(toa), (5,.6), fontsize=6)
        plt.plot(exact, linestyle='steps-post', linewidth=.5)
        plt.plot(psint, linestyle='steps-post', linewidth=.5)
        plt.gca().set_ylim(0,1.2)
        for tick in plt.gca().yaxis.get_major_ticks():
          tick.label1.set_fontsize(6)
        for tick in plt.gca().xaxis.get_major_ticks():
          tick.label1.set_fontsize(6)
        if fig == 1 : fig += 3
 
        table.add_row([
          iord, fct, toa, psint.max() - 1.1, psint.min() - .1, abs(psint - exact).sum(), 
          math.sqrt(((psint - exact)*(psint - exact)).sum()), 
          time / time_ref
        ])

        print 'closing and deleting the file...'
        nc.close()
#        if os.path.exists(file) : os.unlink(file)
print table.draw()
open('tab1.txt','w').write(table.draw() + '\n')
plt.annotate(table.draw(), (.4,.73), xycoords='figure fraction', fontsize=6, family='monospace')
plt.savefig('fig1.svg')
